use super::*;

#[test]
fn checked_order_rejects_values_outside_the_public_range() -> Result<(), &'static str> {
    let Err(error) = ContinuityOrder::new(MAX_CONTINUITY_ORDER + 1) else {
        return Err("orders above G4 must be rejected");
    };

    assert_eq!(error.requested(), 5);
    assert_eq!(error.maximum(), MAX_CONTINUITY_ORDER);
    assert!(!ContinuityOrder::G0.is_experimental());
    assert!(!ContinuityOrder::G3.is_experimental());
    assert!(ContinuityOrder::G4.is_experimental());
    Ok(())
}

#[test]
fn checked_order_conversions_preserve_the_order() {
    for (order, checked) in [
        ContinuityOrder::G0,
        ContinuityOrder::G1,
        ContinuityOrder::G2,
        ContinuityOrder::G3,
        ContinuityOrder::G4,
    ]
    .into_iter()
    .enumerate()
    {
        assert_eq!(ContinuityOrder::try_from(order), Ok(checked));
        assert_eq!(usize::from(checked), order);
        assert_eq!(checked.as_usize(), order);
    }
}

#[test]
fn capability_reports_preserve_every_side_and_typed_support() -> Result<(), &'static str> {
    for side in [
        BoundarySide::MinU,
        BoundarySide::MaxU,
        BoundarySide::MinV,
        BoundarySide::MaxV,
    ] {
        let Ok(supported) = SurfaceContinuityCapability::try_supported_through(
            side,
            ContinuityOrder::G3,
            ContinuityOrder::G4,
        ) else {
            return Err("G4 must be a valid maximum for a G3 request");
        };
        let Ok(unsupported) = SurfaceContinuityCapability::try_unsupported(
            side,
            ContinuityOrder::G4,
            UnsupportedContinuityCapability::InsufficientDegree {
                available: 3,
                required: 4,
            },
            Some(ContinuityOrder::G3),
        ) else {
            return Err("G3 must be a valid maximum for an unsupported G4 request");
        };

        assert_eq!(supported.side(), side);
        assert_eq!(supported.requested(), ContinuityOrder::G3);
        assert_eq!(
            supported.maximum_supported_order(),
            Some(ContinuityOrder::G4)
        );
        assert_eq!(supported.unsupported_reason(), None);
        assert_eq!(unsupported.side(), side);
        assert_eq!(unsupported.requested(), ContinuityOrder::G4);
        assert_eq!(
            unsupported.maximum_supported_order(),
            Some(ContinuityOrder::G3)
        );
        assert_eq!(
            unsupported.unsupported_reason(),
            Some(UnsupportedContinuityCapability::InsufficientDegree {
                available: 3,
                required: 4,
            })
        );
    }
    Ok(())
}

#[test]
fn capability_reports_reject_inconsistent_maximum_orders() -> Result<(), &'static str> {
    let supported = SurfaceContinuityCapability::try_supported_through(
        BoundarySide::MinU,
        ContinuityOrder::G3,
        ContinuityOrder::G2,
    );
    let Err(supported) = supported else {
        return Err("a supported maximum must cover the requested order");
    };
    let unsupported = SurfaceContinuityCapability::try_unsupported(
        BoundarySide::MinU,
        ContinuityOrder::G3,
        UnsupportedContinuityCapability::InsufficientDegree {
            available: 2,
            required: 3,
        },
        Some(ContinuityOrder::G3),
    );
    let Err(unsupported) = unsupported else {
        return Err("an unsupported maximum must be below the requested order");
    };

    assert_eq!(supported.requested(), ContinuityOrder::G3);
    assert_eq!(supported.maximum(), ContinuityOrder::G2);
    assert_eq!(unsupported.requested(), ContinuityOrder::G3);
    assert_eq!(unsupported.maximum(), ContinuityOrder::G3);
    Ok(())
}
