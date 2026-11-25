//! Dense indexed active intervals, preserving swap-remove iteration order.

pub(super) struct ActiveSweep {
    active: Vec<usize>,
    slots: Vec<usize>,
}

impl ActiveSweep {
    const INACTIVE: usize = usize::MAX;

    pub(super) fn new(count: usize) -> Self {
        Self {
            active: Vec::new(),
            slots: vec![Self::INACTIVE; count],
        }
    }

    pub(super) fn insert(&mut self, index: usize) -> bool {
        if self.slots[index] != Self::INACTIVE {
            return false;
        }
        self.slots[index] = self.active.len();
        self.active.push(index);
        true
    }

    pub(super) fn remove(&mut self, index: usize) -> bool {
        let slot = std::mem::replace(&mut self.slots[index], Self::INACTIVE);
        if slot == Self::INACTIVE {
            return false;
        }
        self.active.swap_remove(slot);
        if let Some(&moved) = self.active.get(slot) {
            self.slots[moved] = slot;
        }
        true
    }

    pub(super) fn indices(&self) -> &[usize] { &self.active }
}

#[test]
fn indexed_sweep_preserves_swap_remove_order() {
    let mut sweep = ActiveSweep::new(19);
    let mut reference = Vec::new();
    for round in 0..76 {
        let index = (round * 7) % 19;
        if let Some(slot) = reference.iter().position(|&active| active == index) {
            reference.swap_remove(slot);
            assert!(sweep.remove(index));
            assert!(!sweep.remove(index));
        } else {
            reference.push(index);
            assert!(sweep.insert(index));
            assert!(!sweep.insert(index));
        }
        assert_eq!(sweep.indices(), reference);
    }
    assert!(sweep.indices().is_empty());
}
