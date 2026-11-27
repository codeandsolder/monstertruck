//! Outward-rounded interval arithmetic for complete rational chart spans.
#[derive(Clone, Copy, Debug)]
pub(super) struct Interval {
    pub(super) lo: f64,
    pub(super) hi: f64,
}

impl Interval {
    pub(super) const fn point(x: f64) -> Self { Self { lo: x, hi: x } }
    pub(super) const fn range(a: f64, b: f64) -> Self {
        Self {
            lo: a.min(b),
            hi: a.max(b),
        }
    }
    pub(super) fn finite(self) -> bool {
        self.lo.is_finite() && self.hi.is_finite() && self.lo <= self.hi
    }
    pub(super) const fn union(self, other: Self) -> Self {
        Self {
            lo: self.lo.min(other.lo),
            hi: self.hi.max(other.hi),
        }
    }
    pub(super) fn intersection(self, other: Self) -> Option<Self> {
        let result = Self {
            lo: self.lo.max(other.lo),
            hi: self.hi.min(other.hi),
        };
        result.finite().then_some(result)
    }
    pub(super) fn add(self, other: Self) -> Self {
        Self {
            lo: (self.lo + other.lo).next_down(),
            hi: (self.hi + other.hi).next_up(),
        }
    }
    pub(super) fn sub(self, other: Self) -> Self {
        Self {
            lo: (self.lo - other.hi).next_down(),
            hi: (self.hi - other.lo).next_up(),
        }
    }
    pub(super) fn mul(self, other: Self) -> Self {
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
    pub(super) fn scale(self, x: f64) -> Self { self.mul(Self::point(x)) }
    pub(super) fn div(self, other: Self) -> Option<Self> {
        if !self.finite() || !other.finite() || other.lo <= 0.0 && other.hi >= 0.0 {
            return None;
        }
        let result = self.mul(Self {
            lo: (1.0 / other.hi).next_down(),
            hi: (1.0 / other.lo).next_up(),
        });
        result.finite().then_some(result)
    }
    pub(super) fn square(self) -> Self {
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
    pub(super) fn sqrt(self) -> Option<Self> {
        if !self.finite() || self.hi < 0.0 {
            return None;
        }
        Some(Self {
            lo: self.lo.max(0.0).sqrt().next_down().max(0.0),
            hi: self.hi.sqrt().next_up(),
        })
    }
    pub(super) const fn mid(self) -> f64 { self.lo.midpoint(self.hi) }
    pub(super) const fn magnitude(self) -> f64 { self.lo.abs().max(self.hi.abs()) }
    pub(super) fn width(self) -> f64 { (self.hi - self.lo).next_up() }
    pub(super) fn directed(self) -> bool { self.lo > 0.0 || self.hi < 0.0 }
    #[cfg(test)]
    pub(super) fn contains(self, x: f64) -> bool { self.lo <= x && x <= self.hi }
}

// A narrow rectangle about a nonzero vector has a unique angular lift.
// Bound the angular uncertainty using the Euclidean perturbation radius;
// asin(r) <= 2r for 0<=r<=1/4. The extra angular margin covers libm and
// the floating representation of the chart constants.
pub(super) fn angle(x: Interval, y: Interval) -> Option<Interval> {
    let (a, b) = (x.mid(), y.mid());
    let radius = x
        .sub(Interval::point(a))
        .square()
        .add(y.sub(Interval::point(b)).square())
        .sqrt()?
        .hi;
    let norm = Interval::point(a)
        .square()
        .add(Interval::point(b).square())
        .sqrt()?
        .lo;
    let relative = Interval::point(radius).div(Interval::point(norm))?.hi;
    if relative > 0.25 {
        return None;
    }
    let center = b.atan2(a);
    let error = Interval::point(relative)
        .scale(2.0)
        .add(Interval::point(128.0 * f64::EPSILON))
        .hi;
    Some(Interval::point(center).add(Interval::range(-error, error)))
}
