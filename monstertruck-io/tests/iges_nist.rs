//! End-to-end checks for a small public IGES 4.0 source document.
//!
//! The 28 MB upstream assembly fixture lives in the external regression corpus;
//! this test exercises the decoder, loss reporting and named body refusals
//! without embedding a large binary-sized text asset in this repository.

#![cfg(feature = "iges")]

use monstertruck_io::{Error, iges};

fn fixture() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/iges/nist-ex1.iges")
}

#[test]
fn iges_4_0_reaches_conversion_and_refuses_only_unsupported_wire() {
    let recovered = iges::from_path(fixture()).expect("cadmpeg 0.6 reads IGES 4.0");
    assert!(recovered.bodies.is_empty());
    assert_eq!(recovered.refused.len(), 1, "{}", recovered.summary());
    let refusal = &recovered.refused[0];
    assert_eq!(refusal.kind, "wire");
    assert!(
        matches!(
            refusal.reason,
            Error::UnsupportedBodyKind { kind: "wire", .. }
        ),
        "unexpected refusal: {}",
        refusal.reason
    );
}

#[test]
fn successful_iges_decode_preserves_loss_details() {
    let recovered = iges::from_path(fixture()).expect("IGES 4.0 decodes");
    assert!(
        !recovered.losses.is_empty(),
        "the source deliberately exercises partial decoder losses"
    );
    for loss in &recovered.losses {
        assert!(loss.code.starts_with("iges/"), "{}", loss.code);
        assert!(!loss.message.is_empty());
    }
    assert!(
        recovered
            .losses
            .iter()
            .any(|loss| loss.code == "iges/entity.not-projected"),
        "the codec must report entities omitted during projection"
    );
}

#[test]
fn strict_iges_import_still_names_the_wire_refusal() {
    let err = iges::all_bodies_from_path(fixture())
        .expect_err("a strict import must reject an unsupported wire");
    assert!(
        matches!(err, Error::UnsupportedBodyKind { kind: "wire", .. }),
        "unexpected strict refusal: {err}"
    );
}
