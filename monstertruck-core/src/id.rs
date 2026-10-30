use std::fmt::Debug;
use std::hash::{Hash, Hasher};
use std::marker::PhantomData;

/// Id structure with `Copy`, `Hash` and `Eq` using raw pointers
pub struct Id<T>(usize, PhantomData<T>);

impl<T> Id<T> {
    /// Creates the `Id` from a raw pointer.
    #[must_use]
    #[inline]
    pub fn new(ptr: *const T) -> Self { Self(ptr as usize, PhantomData) }
}

impl<T> Clone for Id<T> {
    #[inline]
    fn clone(&self) -> Self { *self }
}

impl<T> Copy for Id<T> {}

impl<T> Hash for Id<T> {
    #[inline]
    fn hash<H: Hasher>(&self, state: &mut H) { self.0.hash(state); }
}

impl<T> PartialEq for Id<T> {
    #[inline]
    fn eq(&self, other: &Self) -> bool { self.0 == other.0 }
}

impl<T> Eq for Id<T> {}

impl<T> Debug for Id<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        f.write_fmt(format_args!("0x{:x}", self.0))
    }
}

/// Persistent topology element identifier.
///
/// Unlike [`Id`] (derived from Arc pointer addresses, changes on
/// deserialization), `StableId` survives serialization round-trips.
/// Values are unique within a Solid but not globally unique.
#[derive(
    Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, serde::Serialize, serde::Deserialize,
)]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
#[cfg_attr(feature = "rkyv", rkyv(derive(PartialEq, Eq, PartialOrd, Ord, Hash)))]
pub struct StableId(u64);

impl StableId {
    /// Reserved sentinel for unassigned IDs.
    pub const UNASSIGNED: Self = Self(0);

    /// Create a new `StableId`. Zero denotes [`Self::UNASSIGNED`].
    #[must_use]
    pub const fn new(val: u64) -> Self { Self(val) }

    /// Get the raw u64 value.
    #[must_use]
    pub const fn raw(self) -> u64 { self.0 }

    /// Whether this ID has been assigned.
    #[must_use]
    pub const fn is_assigned(self) -> bool { self.0 != 0 }
}

impl Default for StableId {
    fn default() -> Self { Self::UNASSIGNED }
}

impl std::fmt::Display for StableId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result { write!(f, "#{}", self.0) }
}

/// Error returned when a [`StableIdAllocator`] has no identifiers left.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct StableIdExhausted;

impl std::fmt::Display for StableIdExhausted {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("stable ID space is exhausted")
    }
}

impl std::error::Error for StableIdExhausted {}

/// Monotonic counter for assigning fresh [`StableId`]s within a Solid.
#[derive(Clone, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
pub struct StableIdAllocator {
    next: u64,
}

impl StableIdAllocator {
    /// Create a new allocator starting at 1.
    #[must_use]
    pub const fn new() -> Self { Self { next: 1 } }

    /// Allocates the next [`StableId`].
    ///
    /// # Errors
    /// Returns [`StableIdExhausted`] after the `u64` identifier space has been consumed.
    pub fn allocate(&mut self) -> Result<StableId, StableIdExhausted> {
        if self.next == 0 {
            return Err(StableIdExhausted);
        }
        let id = StableId::new(self.next);
        self.next = self.next.checked_add(1).unwrap_or(0);
        Ok(id)
    }

    /// Advances the counter past `ceil` to avoid collisions.
    ///
    /// # Errors
    /// Returns [`StableIdExhausted`] when no identifier exists above `ceil`.
    pub fn ensure_above(&mut self, ceil: u64) -> Result<(), StableIdExhausted> {
        if self.next == 0 {
            return Err(StableIdExhausted);
        }
        if self.next <= ceil {
            self.next = ceil.checked_add(1).ok_or(StableIdExhausted)?;
        }
        Ok(())
    }

    /// Verifies that `count` fresh identifiers can be allocated without mutating the allocator.
    ///
    /// # Errors
    /// Returns [`StableIdExhausted`] when fewer than `count` identifiers remain.
    pub fn ensure_capacity(&self, count: usize) -> Result<(), StableIdExhausted> {
        if count == 0 {
            return Ok(());
        }
        if self.next == 0 {
            return Err(StableIdExhausted);
        }
        let required = u64::try_from(count).map_err(|_| StableIdExhausted)?;
        let available = u64::MAX - self.next + 1;
        (required <= available)
            .then_some(())
            .ok_or(StableIdExhausted)
    }

    /// Current counter value (next ID that will be allocated).
    #[must_use]
    pub const fn peek(&self) -> u64 { self.next }
}

impl Default for StableIdAllocator {
    fn default() -> Self { Self::new() }
}

#[test]
fn debug_formats_id_as_pointer_address() {
    let x: f64 = 3.0;
    let id = Id::new(&raw const x);
    let a = format!("{id:?}");
    let b = format!("{:p}", &x);
    assert_eq!(a, b);
}

#[test]
fn stable_id_basics() {
    let id = StableId::new(42);
    assert_eq!(id.raw(), 42);
    assert!(id.is_assigned());
    assert!(!StableId::UNASSIGNED.is_assigned());
}

#[test]
fn stable_id_allocator() -> Result<(), StableIdExhausted> {
    let mut alloc = StableIdAllocator::new();
    let a = alloc.allocate()?;
    let b = alloc.allocate()?;
    assert_ne!(a, b);
    assert_eq!(a.raw(), 1);
    assert_eq!(b.raw(), 2);
    Ok(())
}

#[test]
fn stable_id_allocator_exhaustion_is_explicit() -> Result<(), StableIdExhausted> {
    assert_eq!(StableId::new(0), StableId::UNASSIGNED);
    let mut alloc = StableIdAllocator::new();
    alloc.ensure_above(u64::MAX - 1)?;
    alloc.ensure_capacity(1)?;
    assert_eq!(alloc.ensure_capacity(2), Err(StableIdExhausted));
    assert_eq!(alloc.allocate()?.raw(), u64::MAX);
    assert_eq!(alloc.allocate(), Err(StableIdExhausted));
    assert_eq!(alloc.peek(), 0);
    Ok(())
}

#[test]
fn failed_ensure_above_does_not_poison_allocator() {
    let mut alloc = StableIdAllocator::new();
    assert_eq!(alloc.ensure_above(u64::MAX), Err(StableIdExhausted));
    assert_eq!(alloc.peek(), 1);
}
