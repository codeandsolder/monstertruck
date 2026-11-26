use anyhow::{Context, Result, bail, ensure};
use flate2::read::GzDecoder;
use futures_util::StreamExt as _;
use serde::Deserialize;
use sha2::{Digest, Sha256};
use std::{
    fs::{self, File},
    io::{self, BufRead, BufReader, Read, Write},
    path::{Component, Path, PathBuf},
    time::Duration,
};
use tempfile::NamedTempFile;
use zip::ZipArchive;

#[derive(Debug, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Archive {
    Zip,
}

#[derive(Debug, Deserialize)]
pub struct Fixture {
    pub filename: Option<String>,
    sha256: String,
    url: String,
    #[serde(default)]
    fallback_urls: Vec<String>,
    bundled_gzip: Option<PathBuf>,
    archive: Option<Archive>,
    archive_filename: Option<String>,
    archive_sha256: Option<String>,
    bundled_archive: Option<PathBuf>,
    member: Option<String>,
}

fn digest(path: &Path) -> Result<String> {
    let mut file = BufReader::new(File::open(path)?);
    let mut hasher = Sha256::new();
    loop {
        let bytes = file.fill_buf()?;
        if bytes.is_empty() {
            break;
        }
        hasher.update(bytes);
        let length = bytes.len();
        file.consume(length);
    }
    Ok(hex::encode(hasher.finalize()))
}

fn validate_digest(expected: &str) -> Result<()> {
    ensure!(
        expected.len() == 64 && expected.bytes().all(|byte| byte.is_ascii_hexdigit()),
        "invalid SHA-256 digest {expected:?}"
    );
    Ok(())
}

fn verify(path: &Path, expected: &str) -> Result<()> {
    validate_digest(expected)?;
    let actual = digest(path)?;
    ensure!(
        actual.eq_ignore_ascii_case(expected),
        "{}: SHA-256 mismatch: expected {expected}, got {actual}",
        path.display()
    );
    Ok(())
}

fn cached(path: &Path, expected: Option<&str>) -> Result<bool> {
    if let Some(expected) = expected {
        validate_digest(expected)?;
    }
    if !path.try_exists()? {
        return Ok(false);
    }
    ensure!(path.is_file(), "{} is not a file", path.display());
    expected.map_or(Ok(true), |expected| {
        Ok(digest(path)?.eq_ignore_ascii_case(expected))
    })
}

pub fn safe_filename(name: &str) -> Result<&str> {
    let mut components = Path::new(name).components();
    ensure!(
        !name.is_empty()
            && !name.contains(['/', '\\', ':'])
            && matches!(components.next(), Some(Component::Normal(_)))
            && components.next().is_none(),
        "unsafe fixture or case filename {name:?}"
    );
    Ok(name)
}

fn bundled_path(base: &Path, relative: &Path) -> Result<PathBuf> {
    ensure!(
        !relative.as_os_str().is_empty()
            && relative
                .components()
                .all(|part| matches!(part, Component::Normal(_))),
        "bundled fixture path must be relative to the manifest directory"
    );
    Ok(base.join(relative))
}

fn install_verified(
    destination: &Path,
    mut source: impl Read,
    expected: Option<&str>,
) -> Result<()> {
    let directory = destination
        .parent()
        .context("destination has no directory")?;
    let mut temporary = NamedTempFile::new_in(directory)?;
    io::copy(&mut source, &mut temporary)?;
    if let Some(expected) = expected {
        verify(temporary.path(), expected)?;
    }
    temporary.persist(destination)?;
    Ok(())
}

async fn download_one(
    client: &wreq::Client,
    url: &str,
    destination: &Path,
    expected: Option<&str>,
) -> Result<()> {
    let response = client
        .get(url)
        .header("User-Agent", "step-redox-validation/0.1")
        .send()
        .await?
        .error_for_status()?;
    let directory = destination
        .parent()
        .context("destination has no directory")?;
    let mut temporary = NamedTempFile::new_in(directory)?;
    let mut stream = std::pin::pin!(response.bytes_stream());
    while let Some(chunk) = stream.next().await {
        temporary.write_all(&chunk?)?;
    }
    if let Some(expected) = expected {
        verify(temporary.path(), expected)?;
    }
    temporary.persist(destination)?;
    Ok(())
}

fn download(spec: &Fixture, destination: &Path, expected: Option<&str>) -> Result<()> {
    // Only the network acquisition boundary is asynchronous; the validator's
    // CPU-intensive geometry and cached-file paths remain synchronous.
    let runtime = tokio::runtime::Builder::new_current_thread()
        .enable_all()
        .build()?;
    runtime.block_on(async {
        let client = wreq::Client::builder()
            .timeout(Duration::from_secs(120))
            .redirect(wreq::redirect::Policy::limited(10))
            .build()?;
        let mut errors = Vec::new();
        for url in std::iter::once(&spec.url).chain(&spec.fallback_urls) {
            match download_one(&client, url, destination, expected).await {
                Ok(()) => return Ok(()),
                Err(error) => errors.push(format!("{url}: {error:#}")),
            }
        }
        bail!("all fixture download URLs failed: {}", errors.join(" | "))
    })
}

pub fn fetch(name: &str, spec: &Fixture, base: &Path, cache: &Path) -> Result<PathBuf> {
    fs::create_dir_all(cache)?;
    let filename = spec
        .filename
        .clone()
        .unwrap_or_else(|| format!("{name}.step"));
    let path = cache.join(safe_filename(&filename)?);
    if cached(&path, Some(&spec.sha256))? {
        return Ok(path);
    }

    if let Some(relative) = &spec.bundled_gzip {
        let bundled = bundled_path(base, relative)?;
        if bundled.try_exists()? {
            install_verified(
                &path,
                GzDecoder::new(File::open(&bundled)?),
                Some(&spec.sha256),
            )
            .with_context(|| format!("extract {}", bundled.display()))?;
            return Ok(path);
        }
    }

    match spec.archive {
        Some(Archive::Zip) => {
            let filename = spec
                .archive_filename
                .clone()
                .unwrap_or_else(|| format!("{filename}.zip"));
            let archive = cache.join(safe_filename(&filename)?);
            ensure!(
                archive != path,
                "fixture archive and STEP filenames coincide"
            );
            let expected = spec.archive_sha256.as_deref();
            if !cached(&archive, expected)? {
                let bundled = spec
                    .bundled_archive
                    .as_deref()
                    .map(|relative| bundled_path(base, relative))
                    .transpose()?;
                if let Some(bundled) = bundled.filter(|path| path.exists()) {
                    install_verified(&archive, File::open(bundled)?, expected)?;
                } else {
                    download(spec, &archive, expected)?;
                }
            }
            let member = spec
                .member
                .as_deref()
                .context("ZIP fixture has no member")?;
            let mut archive = ZipArchive::new(File::open(archive)?)?;
            let entry = archive.by_name(member)?;
            ensure!(entry.is_file(), "ZIP fixture member must be a file");
            // Read exactly the named member; never extract archive paths onto disk.
            install_verified(&path, entry, Some(&spec.sha256))?;
        }
        None => download(spec, &path, Some(&spec.sha256))?,
    }
    Ok(path)
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::{Compression, write::GzEncoder};

    #[test]
    fn verifies_and_repairs_cache_without_partial_publication() -> Result<()> {
        let root = tempfile::tempdir()?;
        let cache = root.path().join("cache");
        let bytes = b"ISO-10303-21; END-ISO-10303-21;";
        let expected = hex::encode(Sha256::digest(bytes));
        let gzip = root.path().join("part.step.gz");
        let mut encoder = GzEncoder::new(File::create(&gzip)?, Compression::default());
        encoder.write_all(bytes)?;
        encoder.finish()?;
        let spec: Fixture = serde_json::from_value(serde_json::json!({
            "url": "https://invalid.example", "filename": "part.step",
            "sha256": expected, "bundled_gzip": "part.step.gz"
        }))?;
        let path = fetch("part", &spec, root.path(), &cache)?;
        assert_eq!(fs::read(&path)?, bytes);
        fs::write(&path, b"old corrupt bytes")?;
        install_verified(&path, &b"bad download"[..], Some(&spec.sha256)).expect_err("bad digest");
        assert_eq!(fs::read(&path)?, b"old corrupt bytes");
        fetch("part", &spec, root.path(), &cache)?;
        assert_eq!(fs::read(path)?, bytes);
        assert_eq!(fs::read_dir(cache)?.count(), 1);
        Ok(())
    }

    #[test]
    fn rejects_escaping_paths() {
        for name in ["", ".", "..", "../x", "a/b", "a\\b", "C:drive"] {
            assert!(safe_filename(name).is_err(), "{name}");
        }
        assert!(bundled_path(Path::new("root"), Path::new("../x")).is_err());
        assert!(bundled_path(Path::new("root"), Path::new("/absolute")).is_err());
    }
}
