//! Outward-rounded interval arithmetic for complete rational chart spans.
/// Closed interval enclosing a real scalar. Arithmetic rounds outward.
#[derive(Clone, Copy, Debug)]
#[must_use = "interval operations return a new enclosure"]
pub struct Interval {
    /// Inclusive lower endpoint.
    pub lo: f64,
    /// Inclusive upper endpoint.
    pub hi: f64,
}

impl Interval {
    /// Enclose an exactly represented floating-point number.
    #[must_use = "inspect the computed enclosure or predicate"]
    pub const fn point(x: f64) -> Self { Self { lo: x, hi: x } }
    /// Enclose the interval between two endpoints.
    #[must_use = "inspect the computed enclosure or predicate"]
    pub const fn range(a: f64, b: f64) -> Self {
        Self {
            lo: a.min(b),
            hi: a.max(b),
        }
    }
    /// Check for a finite nonempty enclosure.
    #[must_use = "inspect the computed enclosure or predicate"]
    pub const fn finite(self) -> bool {
        self.lo.is_finite() && self.hi.is_finite() && self.lo <= self.hi
    }
    /// Smallest closed interval containing both enclosures.
    #[must_use = "inspect the computed enclosure or predicate"]
    pub const fn union(self, other: Self) -> Self {
        Self {
            lo: self.lo.min(other.lo),
            hi: self.hi.max(other.hi),
        }
    }
    /// Intersect two enclosures; reject an empty or nonfinite result.
    #[must_use = "inspect the computed enclosure or predicate"]
    pub fn intersection(self, other: Self) -> Option<Self> {
        let result = Self {
            lo: self.lo.max(other.lo),
            hi: self.hi.min(other.hi),
        };
        result.finite().then_some(result)
    }
    /// Outward-rounded sum.
    #[must_use = "inspect the computed enclosure or predicate"]
    pub fn plus(self, other: Self) -> Self {
        Self {
            lo: (self.lo + other.lo).next_down(),
            hi: (self.hi + other.hi).next_up(),
        }
    }
    /// Outward-rounded difference.
    #[must_use = "inspect the computed enclosure or predicate"]
    pub fn minus(self, other: Self) -> Self {
        Self {
            lo: (self.lo - other.hi).next_down(),
            hi: (self.hi - other.lo).next_up(),
        }
    }
    /// Outward-rounded product.
    #[must_use = "inspect the computed enclosure or predicate"]
    pub fn times(self, other: Self) -> Self {
        let products = [
            self.lo * other.lo,
            self.lo * other.hi,
            self.hi * other.lo,
            self.hi * other.hi,
        ];
        if products.iter().any(|x| !x.is_finite()) {
            return Self {
                lo: f64::NEG_INFINITY,
                hi: f64::INFINITY,
            };
        }
        Self {
            lo: products
                .into_iter()
                .fold(f64::INFINITY, f64::min)
                .next_down(),
            hi: products
                .into_iter()
                .fold(f64::NEG_INFINITY, f64::max)
                .next_up(),
        }
    }
    /// Multiply by an exactly represented scalar.
    #[must_use = "inspect the computed enclosure or predicate"]
    pub fn scale(self, x: f64) -> Self { self.times(Self::point(x)) }
    /// Outward-rounded quotient, refusing a divisor containing zero.
    #[must_use = "inspect the computed enclosure or predicate"]
    pub fn checked_div(self, other: Self) -> Option<Self> {
        if !self.finite() || !other.finite() || other.lo <= 0.0 && other.hi >= 0.0 {
            return None;
        }
        let result = self.times(Self {
            lo: (1.0 / other.hi).next_down(),
            hi: (1.0 / other.lo).next_up(),
        });
        result.finite().then_some(result)
    }
    /// Outward-rounded square.
    #[must_use = "inspect the computed enclosure or predicate"]
    pub fn square(self) -> Self {
        let [a, b] = [self.lo * self.lo, self.hi * self.hi];
        let lo = if self.lo <= 0.0 && self.hi >= 0.0 {
            0.0
        } else {
            a.min(b).next_down().max(0.0)
        };
        Self {
            lo,
            hi: a.max(b).next_up(),
        }
    }
    /// Outward-rounded square root over the nonnegative part.
    #[must_use = "inspect the computed enclosure or predicate"]
    pub fn sqrt(self) -> Option<Self> {
        if !self.finite() || self.hi < 0.0 {
            return None;
        }
        Some(Self {
            lo: self.lo.max(0.0).sqrt().next_down().max(0.0),
            hi: self.hi.sqrt().next_up(),
        })
    }
    /// Floating-point midpoint used to choose refinement witnesses.
    #[must_use = "inspect the computed enclosure or predicate"]
    pub const fn mid(self) -> f64 { self.lo.midpoint(self.hi) }
    /// Upper bound on absolute value.
    #[must_use = "inspect the computed enclosure or predicate"]
    pub const fn magnitude(self) -> f64 { self.lo.abs().max(self.hi.abs()) }
    /// Upper bound on the interval width.
    #[must_use = "inspect the computed enclosure or predicate"]
    pub fn width(self) -> f64 { (self.hi - self.lo).next_up() }
    /// Check whether zero is excluded.
    #[must_use = "inspect the computed enclosure or predicate"]
    pub const fn directed(self) -> bool { self.lo > 0.0 || self.hi < 0.0 }
    /// Check whether a represented scalar is enclosed.
    #[must_use = "inspect the computed enclosure or predicate"]
    pub const fn contains(self, x: f64) -> bool { self.lo <= x && x <= self.hi }
}

// A narrow rectangle about a nonzero vector has a unique angular lift.
// Bound the angular uncertainty using the Euclidean perturbation radius;
// asin(r) <= 2r for 0<=r<=1/4. The extra angular margin covers libm and
// the floating representation of the chart constants.
/// Enclose the angle of a nonzero vector in a sufficiently narrow rectangle.
#[must_use = "an unresolved angular lift is not a valid enclosure"]
pub fn angular_enclosure(x: Interval, y: Interval) -> Option<Interval> {
    let (a, b) = (x.mid(), y.mid());
    let radius = x
        .minus(Interval::point(a))
        .square()
        .plus(y.minus(Interval::point(b)).square())
        .sqrt()?
        .hi;
    let norm = Interval::point(a)
        .square()
        .plus(Interval::point(b).square())
        .sqrt()?
        .lo;
    let relative = Interval::point(radius)
        .checked_div(Interval::point(norm))?
        .hi;
    if relative > 0.25 {
        return None;
    }
    let center = b.atan2(a);
    let error = Interval::point(relative)
        .scale(2.0)
        .plus(Interval::point(128.0 * f64::EPSILON))
        .hi;
    Some(Interval::point(center).plus(Interval::range(-error, error)))
}
