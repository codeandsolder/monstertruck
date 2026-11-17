use anyhow::{Context, Result, bail};
use std::path::PathBuf;

/// Shared positional contract for the extrusion-core recovery diagnostics.
pub(super) struct ExtrusionCoreArgs {
    input: PathBuf,
    pub(super) solid_id: u64,
    pub(super) caps: [u64; 2],
    pub(super) output: PathBuf,
}

impl ExtrusionCoreArgs {
    pub(super) fn parse(usage: &str) -> Result<Self> {
        Self::parse_from(std::env::args().skip(1), usage)
    }

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

    pub(super) fn read_input(&self) -> Result<Vec<u8>> {
        std::fs::read(&self.input).with_context(|| format!("read {}", self.input.display()))
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
    fn rejects_missing_excess_and_non_numeric_arguments() {
        assert!(args(&["input.step", "1", "2", "3"]).is_err());
        assert!(args(&["input.step", "1", "2", "3", "out", "extra"]).is_err());
        assert!(args(&["input.step", "nope", "2", "3", "out"]).is_err());
        assert!(args(&["input.step", "1", "-2", "3", "out"]).is_err());
    }
}
