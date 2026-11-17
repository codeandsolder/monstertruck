use anyhow::{Context, Result, bail};
use std::path::PathBuf;

/// Shared positional contract for the extrusion-core recovery diagnostics.
pub struct ExtrusionCoreArgs {
    input: PathBuf,
    pub solid_id: u64,
    pub caps: [u64; 2],
    pub output: PathBuf,
}

impl ExtrusionCoreArgs {
    pub fn parse(usage: &str) -> Result<Self> { Self::parse_from(std::env::args().skip(1), usage) }

    fn parse_from(args: impl IntoIterator<Item = String>, usage: &str) -> Result<Self> {
        let mut args = args.into_iter();
        let input = PathBuf::from(args.next().context(usage.to_owned())?);
        let solid_id = args
            .next()
            .context(usage.to_owned())?
            .parse::<u64>()
            .context("parse SOLID_ID")?;
        let cap_a = args
            .next()
            .context(usage.to_owned())?
            .parse::<u64>()
            .context("parse CAP_A")?;
        let cap_b = args
            .next()
            .context(usage.to_owned())?
            .parse::<u64>()
            .context("parse CAP_B")?;
        let output = PathBuf::from(args.next().context(usage.to_owned())?);
        if args.next().is_some() {
            bail!("{usage}");
        }
        Ok(Self {
            input,
            solid_id,
            caps: [cap_a, cap_b],
            output,
        })
    }

    pub fn read_input(&self) -> Result<Vec<u8>> {
        std::fs::read(&self.input).with_context(|| format!("read {}", self.input.display()))
    }

    /// Resolve an extrusion-core diagnostic using its own recovery and admission rule.
    pub fn recover<T>(
        &self,
        bytes: &[u8],
        description: &str,
        recover: impl FnOnce(&[u8], u64, [u64; 2]) -> Result<Option<T>>,
    ) -> Result<T> {
        recover(bytes, self.solid_id, self.caps)?.with_context(|| {
            let [cap_a, cap_b] = self.caps;
            format!(
                "no {description} for solid #{} caps #{cap_a}/#{cap_b}",
                self.solid_id
            )
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(values: &[&str]) -> Result<ExtrusionCoreArgs> {
        ExtrusionCoreArgs::parse_from(values.iter().map(|s| (*s).to_owned()), "usage: test")
    }

    #[test]
    fn parses_complete_positional_contract() {
        let parsed = args(&["input.step", "42", "12", "34", "out.step"])
            .expect("valid extrusion-core arguments");
        assert_eq!(parsed.input, PathBuf::from("input.step"));
        assert_eq!(parsed.solid_id, 42);
        assert_eq!(parsed.caps, [12, 34]);
        assert_eq!(parsed.output, PathBuf::from("out.step"));
    }

    #[test]
    fn recovery_forwards_source_identity_and_preserves_missing_candidate_diagnostic() {
        let parsed = args(&["source.step", "42", "12", "34", "out.step"]).expect("valid args");
        let result = parsed
            .recover(b"step", "candidate", |bytes, id, caps| {
                assert_eq!(bytes, b"step");
                assert_eq!(id, 42);
                assert_eq!(caps, [12, 34]);
                Ok(Some(7_u8))
            })
            .expect("candidate");
        assert_eq!(result, 7);
        let error = parsed
            .recover(b"", "candidate", |_, _, _| Ok(None::<u8>))
            .expect_err("no candidate");
        assert_eq!(error.to_string(), "no candidate for solid #42 caps #12/#34");
    }

    #[test]
    fn rejects_missing_excess_and_non_numeric_arguments() {
        assert!(args(&["input.step", "1", "2", "3"]).is_err());
        assert!(args(&["input.step", "1", "2", "3", "out", "extra"]).is_err());
        assert!(args(&["input.step", "nope", "2", "3", "out"]).is_err());
        assert!(args(&["input.step", "1", "-2", "3", "out"]).is_err());
    }
}
