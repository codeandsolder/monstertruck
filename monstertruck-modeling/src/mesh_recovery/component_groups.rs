//! Stable disjoint-set grouping for source and candidate incidence.

pub(super) struct Groups {
    parents: Vec<usize>,
}

impl Groups {
    pub(super) fn new(count: usize) -> Self {
        Self {
            parents: (0..count).collect(),
        }
    }

    pub(super) fn root(&mut self, mut id: usize) -> usize {
        while self.parents[id] != id {
            self.parents[id] = self.parents[self.parents[id]];
            id = self.parents[id];
        }
        id
    }

    pub(super) fn join(&mut self, a: usize, b: usize) {
        let a = self.root(a);
        let b = self.root(b);
        self.parents[a.max(b)] = a.min(b);
    }
}
