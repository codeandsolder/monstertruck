//! Durable, version-bound checkpoints for long-running recovery surveys.
use super::{Case, SCHEMA};
use anyhow::{Context, Result, ensure};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::{
    fs::{File, OpenOptions},
    io::{BufRead, BufReader, Read, Seek, SeekFrom, Write},
    path::{Path, PathBuf},
};

#[derive(Debug, PartialEq, Eq, Serialize, Deserialize)]
pub(super) struct Header {
    schema: String,
    pub(super) executable_sha256: String,
    inputs: Vec<PathBuf>,
}

impl Header {
    pub(super) fn current(inputs: &[PathBuf]) -> Result<Self> {
        let executable = std::env::current_exe().context("locate survey executable")?;
        Ok(Self {
            schema: SCHEMA.to_owned(),
            executable_sha256: hash_file(&executable)?.1,
            inputs: inputs.to_vec(),
        })
    }
}

fn hash_file(path: &Path) -> Result<(usize, String)> {
    let mut file = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let mut hasher = Sha256::new();
    let mut buffer = [0_u8; 8192];
    let mut size = 0_usize;
    loop {
        let count = file.read(&mut buffer)?;
        if count == 0 {
            break;
        }
        size = size.checked_add(count).context("input length overflow")?;
        hasher.update(&buffer[..count]);
    }
    Ok((size, hex::encode(hasher.finalize())))
}

fn validate_case(entry: &Case, expected: &Path) -> Result<()> {
    ensure!(
        expected.to_str() == Some(entry.input.as_str()),
        "checkpoint input order differs at {}",
        expected.display()
    );
    match hash_file(expected) {
        Ok((size, sha256)) => ensure!(
            entry.input_bytes == Some(size) && entry.sha256.as_deref() == Some(sha256.as_str()),
            "checkpoint input contents changed: {}",
            expected.display()
        ),
        Err(error) => {
            ensure!(
                entry.input_bytes.is_none()
                    && entry.sha256.is_none()
                    && entry.outcome == super::Outcome::InputError,
                "cannot reverify completed checkpoint input {}: {error:#}",
                expected.display()
            );
        }
    }
    Ok(())
}

fn write_record(file: &mut File, record: &impl Serialize) -> Result<()> {
    let mut bytes = serde_json::to_vec(record)?;
    bytes.push(b'\n');
    file.write_all(&bytes)?;
    // Unlike flush(), this commits each complete result to the filesystem
    // before the next potentially expensive model begins.
    file.sync_data()?;
    Ok(())
}

pub(super) struct Journal {
    file: File,
}

impl Journal {
    pub(super) fn open(path: &Path, header: &Header, resume: bool) -> Result<(Self, Vec<Case>)> {
        let mut file = OpenOptions::new()
            .read(true)
            .write(true)
            .create_new(!resume)
            .open(path)
            .with_context(|| {
                format!(
                    "{} checkpoint {} (use a fresh output or explicit --resume)",
                    if resume { "open" } else { "create" },
                    path.display()
                )
            })?;
        file.try_lock()
            .with_context(|| format!("lock checkpoint {}", path.display()))?;
        let cases = if resume {
            Self::read_prefix(&mut file, header)?
        } else {
            write_record(&mut file, header)?;
            Vec::new()
        };
        Ok((Self { file }, cases))
    }

    fn read_prefix(file: &mut File, expected: &Header) -> Result<Vec<Case>> {
        let mut reader = BufReader::new(&mut *file);
        let mut line = Vec::new();
        let header_length = reader.read_until(b'\n', &mut line)?;
        ensure!(
            line.ends_with(b"\n"),
            "checkpoint header is missing or incomplete"
        );
        let stored: Header = serde_json::from_slice(&line).context("parse checkpoint header")?;
        ensure!(
            &stored == expected,
            "checkpoint schema, executable digest or ordered input list differs"
        );
        let mut offset = u64::try_from(header_length)?;
        let mut cases = Vec::new();
        loop {
            line.clear();
            let count = reader.read_until(b'\n', &mut line)?;
            if count == 0 || !line.ends_with(b"\n") {
                // An interrupted final write is not a completed record.
                // Malformed *terminated* records instead fail closed.
                break;
            }
            let entry: Case = serde_json::from_slice(&line)
                .with_context(|| format!("parse checkpoint case {}", cases.len() + 1))?;
            let input = expected
                .inputs
                .get(cases.len())
                .context("checkpoint has more records than requested inputs")?;
            validate_case(&entry, input)?;
            offset += u64::try_from(count)?;
            cases.push(entry);
        }
        drop(reader);
        file.set_len(offset)?;
        file.seek(SeekFrom::Start(offset))?;
        file.sync_data()?;
        Ok(cases)
    }

    pub(super) fn append(&mut self, entry: &Case) -> Result<()> {
        write_record(&mut self.file, entry)
    }
}

pub(super) fn write_report(path: &Path, bytes: &[u8]) -> Result<()> {
    let directory = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let mut file = tempfile::NamedTempFile::new_in(directory)?;
    file.write_all(bytes)?;
    file.as_file().sync_all()?;
    // Publish only a completely written report; never overwrite a source or
    // a pre-existing report, including if another process wins a race.
    file.persist_noclobber(path)
        .with_context(|| format!("publish complete report {}", path.display()))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::survey::{Outcome, recover_case};
    use step_redox::cad_kernel::monstertruck::MonstertruckKernel;

    fn header(path: &Path) -> Header {
        Header {
            schema: SCHEMA.to_owned(),
            executable_sha256: "test-build".to_owned(),
            inputs: vec![path.to_path_buf(), path.with_extension("next.step")],
        }
    }

    #[test]
    fn resumes_exact_prefix_and_truncates_only_torn_final_record() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let input = dir.path().join("source.step");
        std::fs::write(&input, b"not a STEP document")?;
        let checkpoint = dir.path().join("report.cases.jsonl");
        let header = header(&input);
        let (mut journal, cases) = Journal::open(&checkpoint, &header, false)?;
        assert!(cases.is_empty());
        assert!(Journal::open(&checkpoint, &header, true).is_err());
        let case = recover_case(&input, &MonstertruckKernel::project_standard());
        journal.append(&case)?;
        let complete_len = journal.file.metadata()?.len();
        journal.file.write_all(b"{\"input\":")?;
        drop(journal);
        let (journal, cases) = Journal::open(&checkpoint, &header, true)?;
        assert_eq!(cases.len(), 1);
        assert_eq!(cases[0].outcome, Outcome::InputError);
        assert_eq!(cases[0].sha256, case.sha256);
        assert_eq!(journal.file.metadata()?.len(), complete_len);
        drop(journal);
        std::fs::write(&input, b"changed source")?;
        assert!(Journal::open(&checkpoint, &header, true).is_err());
        assert!(Journal::open(&checkpoint, &header, false).is_err());
        Ok(())
    }

    #[test]
    fn rejects_mismatched_build_and_complete_corrupt_record() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let input = dir.path().join("source.step");
        let checkpoint = dir.path().join("report.cases.jsonl");
        let mut header = header(&input);
        let (mut journal, _) = Journal::open(&checkpoint, &header, false)?;
        journal.file.write_all(b"invalid JSON\n")?;
        drop(journal);
        assert!(Journal::open(&checkpoint, &header, true).is_err());
        header.executable_sha256.push('2');
        assert!(Journal::open(&checkpoint, &header, true).is_err());
        Ok(())
    }

    #[test]
    fn report_publication_never_clobbers_existing_data() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let output = dir.path().join("report.json");
        write_report(&output, b"first")?;
        assert!(write_report(&output, b"second").is_err());
        assert_eq!(std::fs::read(&output)?, b"first");
        Ok(())
    }
}
