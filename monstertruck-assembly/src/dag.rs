use rayon::prelude::*;
use std::collections::BTreeSet;

/// Index for access to the node in the graph
#[derive(Clone, Copy, Debug, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct NodeIndex(usize);

/// Edge of graph
#[derive(Clone, Copy, Debug)]
pub struct EdgeData<EE> {
    to: NodeIndex,
    entity: EE,
}

/// Reference of edge
#[derive(Debug)]
pub struct Edge<'a, EE> {
    nodes: (NodeIndex, NodeIndex),
    entity: &'a EE,
}
impl<'a, EE> Clone for Edge<'a, EE> {
    fn clone(&self) -> Self { *self }
}

impl<'a, EE> Copy for Edge<'a, EE> {}

fn edge_closure<EE>(index: NodeIndex) -> impl Fn(&EdgeData<EE>) -> Edge<'_, EE> {
    move |edge| Edge {
        nodes: (index, edge.to),
        entity: &edge.entity,
    }
}

/// Reference of edge
#[derive(Debug)]
pub struct EdgeMut<'a, EE> {
    nodes: (NodeIndex, NodeIndex),
    entity: &'a mut EE,
}

fn edge_mut_closure<EE>(index: NodeIndex) -> impl Fn(&mut EdgeData<EE>) -> EdgeMut<'_, EE> {
    move |edge| EdgeMut {
        nodes: (index, edge.to),
        entity: &mut edge.entity,
    }
}

/// Node of graph
#[derive(Clone, Debug)]
#[doc(hidden)]
pub struct NodeData<NE, EE> {
    edges: Vec<EdgeData<EE>>,
    parents: BTreeSet<NodeIndex>,
    entity: NE,
}

impl<NE, EE> NodeData<NE, EE> {
    #[inline]
    fn new(entity: NE) -> Self {
        NodeData {
            edges: Default::default(),
            parents: Default::default(),
            entity,
        }
    }
}

/// Node, interface for referencing nodes within a graph
/// # Examples
/// ```
/// use monstertruck_assembly::dag::*;
///
/// // create dag
/// let mut dag = Dag::<usize, usize>::new();
///
/// // create nodes with entity 0, 1 and get indices
/// let a = dag.create_nodes(0..2);
///
/// // create some edges
/// dag.create_edge(a[1], a[0], 0);
/// dag.create_edge(a[1], a[0], 1);
///
/// // get reference of node
/// let Some(node) = dag.node(a[1]) else { return };
/// assert_eq!(node.index(), a[1]);
/// assert_eq!(node.entity(), &1);
///
/// // edge information is stored in the root node.
/// assert_eq!(node.edges().len(), 2);
/// let Some(edge) = node.edge(1) else { return };
/// assert_eq!(edge.entity(), &1);
///
/// // one can follow the parent index.
/// assert!(dag.node(a[0]).is_some_and(|node| node.parents().contains(&a[1])));
/// assert!(dag.node(a[1]).is_some_and(|node| node.parents().is_empty()));
/// ```
#[derive(Debug)]
pub struct Node<'a, NE, EE> {
    index: NodeIndex,
    data: &'a NodeData<NE, EE>,
}

impl<'a, NE, EE> Clone for Node<'a, NE, EE> {
    fn clone(&self) -> Self { *self }
}

impl<'a, NE, EE> Copy for Node<'a, NE, EE> {}

fn node<NE, EE>((index, data): (usize, &NodeData<NE, EE>)) -> Node<'_, NE, EE> {
    Node {
        index: NodeIndex(index),
        data,
    }
}

/// Mutable node, interface for mutable referencing nodes within a graph
/// # Examples
/// ```
/// use monstertruck_assembly::dag::*;
///
/// // create dag
/// let mut dag = Dag::<usize, usize>::new();
///
/// // create nodes with entity 0, 1 and get indices
/// let a = dag.create_nodes(0..2);
///
/// // create some edges
/// dag.create_edge(a[1], a[0], 0);
/// dag.create_edge(a[1], a[0], 1);
///
/// // get reference of node
/// if let Some(mut node) = dag.node_mut(a[1]) {
///     *node.entity() = 2;
/// }
///
/// assert!(dag.node(a[1]).is_some_and(|node| node.entity() == &2));
///
/// // edge information is stored in the root node.
/// if let Some(mut node) = dag.node_mut(a[1]) {
///     if let Some(mut edge) = node.edge(1) {
///         *edge.entity() = 3;
///     }
/// }
/// assert!(dag
///     .node(a[1])
///     .and_then(|node| node.edge(1))
///     .is_some_and(|edge| edge.entity() == &3));
/// ```
#[derive(Debug)]
pub struct NodeMut<'a, NE, EE> {
    /// Index of the node
    index: NodeIndex,
    /// Reference of the node
    data: &'a mut NodeData<NE, EE>,
}

fn node_mut<NE, EE>((index, data): (usize, &mut NodeData<NE, EE>)) -> NodeMut<'_, NE, EE> {
    NodeMut {
        index: NodeIndex(index),
        data,
    }
}

/// DAG
/// # Examples
/// ```
/// use monstertruck_assembly::dag::*;
///
/// // create dag
/// let mut dag = Dag::<usize, usize>::new();
///
/// // create nodes with entity 0, 1,.., 9 and get indices
/// let a = dag.create_nodes(0..10);
///
/// // create some edges
/// dag.create_edge(a[3], a[1], 0);
/// dag.create_edge(a[6], a[4], 1);
/// dag.create_edge(a[3], a[2], 2);
///
/// // get reference of node
/// let Some(node3) = dag.node(a[3]) else { return };
/// // the entity of the node with index 3 is 3.
/// assert_eq!(node3.entity(), &3);
/// // the node with index 3 has 2 edges.
/// assert_eq!(node3.edges().len(), 2);
/// ```
#[derive(Clone, Debug)]
pub struct Dag<NE, EE> {
    nodes: Vec<NodeData<NE, EE>>,
}

/// Path of graph
#[derive(Clone, Debug)]
pub struct Path<'a, NE, EE> {
    nodes: Vec<Node<'a, NE, EE>>,
    edges: Vec<Edge<'a, EE>>,
}

impl<NE, EE> Default for Dag<NE, EE> {
    #[inline]
    fn default() -> Self { Self { nodes: Vec::new() } }
}

impl<NE, EE> Dag<NE, EE> {
    /// Constructor
    #[inline]
    pub const fn new() -> Self { Self { nodes: Vec::new() } }

    /// Constructs a new dag with capacity for `n` values pre-allocated.
    #[inline]
    pub fn with_capacity(n: usize) -> Self {
        Self {
            nodes: Vec::with_capacity(n),
        }
    }

    /// Returns `true` if the graph contains no nodes.
    #[inline]
    pub const fn is_empty(&self) -> bool { self.nodes.is_empty() }

    /// Returns the number of the nodes.
    /// # Examples
    /// ```
    /// use monstertruck_assembly::dag::*;
    /// let mut dag = Dag::<(), ()>::new();
    /// let _ = dag.create_nodes([(); 5]);
    /// assert_eq!(dag.len(), 5);
    /// ```
    #[inline]
    pub const fn len(&self) -> usize { self.nodes.len() }

    /// Returns an iterator over all node indices.
    /// # Examples
    /// ```
    /// use monstertruck_assembly::dag::*;
    /// let mut dag = Dag::<(), ()>::new();
    /// let a = dag.create_nodes([(); 5]);
    /// dag.create_edge(a[0], a[1], ());
    /// dag.create_edge(a[0], a[2], ());
    /// dag.create_edge(a[2], a[3], ());
    /// dag.create_edge(a[4], a[3], ());
    /// let nodes = dag.node_indices().collect::<Vec<_>>();
    /// assert_eq!(nodes.len(), 5);
    /// ```
    #[inline]
    pub fn node_indices(&self) -> NodeIndices { (0..self.nodes.len()).map(NodeIndex) }

    /// Returns a parallel iterator over all node indices.
    /// # Examples
    /// ```
    /// use monstertruck_assembly::dag::*;
    /// use rayon::prelude::*;
    /// let mut dag = Dag::<(), ()>::new();
    /// let a = dag.create_nodes([(); 5]);
    /// dag.create_edge(a[0], a[1], ());
    /// dag.create_edge(a[0], a[2], ());
    /// dag.create_edge(a[2], a[3], ());
    /// dag.create_edge(a[4], a[3], ());
    /// let nodes = dag.par_node_indices().collect::<Vec<_>>();
    /// assert_eq!(nodes.len(), 5);
    /// ```
    #[inline]
    pub fn par_node_indices(&self) -> ParNodeIndices {
        (0..self.nodes.len()).into_par_iter().map(NodeIndex)
    }

    /// Returns all nodes.
    /// # Examples
    /// ```
    /// use monstertruck_assembly::dag::*;
    /// let mut dag = Dag::<(), ()>::new();
    /// let a = dag.create_nodes([(); 5]);
    /// dag.create_edge(a[0], a[1], ());
    /// dag.create_edge(a[0], a[2], ());
    /// dag.create_edge(a[2], a[3], ());
    /// dag.create_edge(a[4], a[3], ());
    /// let nodes = dag.all_nodes().collect::<Vec<_>>();
    /// assert_eq!(nodes.len(), 5);
    /// ```
    #[inline]
    pub fn all_nodes(&self) -> AllNodes<'_, NE, EE> { self.nodes.iter().enumerate().map(node) }

    /// Returns all nodes.
    /// # Examples
    /// ```
    /// use monstertruck_assembly::dag::*;
    /// let mut dag = Dag::<(), ()>::new();
    /// let a = dag.create_nodes([(); 5]);
    /// dag.create_edge(a[0], a[1], ());
    /// dag.create_edge(a[0], a[2], ());
    /// dag.create_edge(a[2], a[3], ());
    /// dag.create_edge(a[4], a[3], ());
    /// let nodes = dag.all_nodes_mut().collect::<Vec<_>>();
    /// assert_eq!(nodes.len(), 5);
    /// ```
    #[inline]
    pub fn all_nodes_mut(&mut self) -> AllNodesMut<'_, NE, EE> {
        self.nodes.iter_mut().enumerate().map(node_mut)
    }

    /// Returns all nodes parallel iterator.
    /// # Examples
    /// ```
    /// use monstertruck_assembly::dag::*;
    /// use rayon::prelude::*;
    /// let mut dag = Dag::<(), ()>::new();
    /// let a = dag.create_nodes([(); 5]);
    /// dag.create_edge(a[0], a[1], ());
    /// dag.create_edge(a[0], a[2], ());
    /// dag.create_edge(a[2], a[3], ());
    /// dag.create_edge(a[4], a[3], ());
    /// let nodes = dag.par_all_nodes().collect::<Vec<_>>();
    /// assert_eq!(nodes.len(), 5);
    /// ```
    #[inline]
    pub fn par_all_nodes(&self) -> ParAllNodes<'_, NE, EE>
    where
        NE: Sync,
        EE: Sync, {
        self.nodes.par_iter().enumerate().map(node)
    }

    /// Returns all nodes parallel iterator.
    /// # Examples
    /// ```
    /// use monstertruck_assembly::dag::*;
    /// use rayon::prelude::*;
    /// let mut dag = Dag::<(), ()>::new();
    /// let a = dag.create_nodes([(); 5]);
    /// dag.create_edge(a[0], a[1], ());
    /// dag.create_edge(a[0], a[2], ());
    /// dag.create_edge(a[2], a[3], ());
    /// dag.create_edge(a[4], a[3], ());
    /// let nodes = dag.par_all_nodes_mut().collect::<Vec<_>>();
    /// assert_eq!(nodes.len(), 5);
    /// ```
    #[inline]
    pub fn par_all_nodes_mut(&mut self) -> ParAllNodesMut<'_, NE, EE>
    where
        NE: Send,
        EE: Send, {
        self.nodes.par_iter_mut().enumerate().map(node_mut)
    }

    /// Returns an iterator on all top level nodes.
    /// # Examples
    /// ```
    /// use monstertruck_assembly::dag::*;
    /// let mut dag = Dag::<usize, ()>::new();
    /// let a = dag.create_nodes(0..5);
    /// dag.create_edge(a[0], a[1], ());
    /// dag.create_edge(a[0], a[2], ());
    /// dag.create_edge(a[2], a[3], ());
    /// dag.create_edge(a[4], a[3], ());
    /// let tops = dag.top_nodes().collect::<Vec<_>>();
    /// assert_eq!(tops.len(), 2);
    /// assert_ne!(tops[0].entity(), tops[1].entity());
    /// ```
    #[inline]
    pub fn top_nodes(&self) -> TopNodes<'_, NE, EE> {
        self.nodes
            .iter()
            .enumerate()
            .filter_map(move |(index, data)| match data.parents.is_empty() {
                true => Some(node((index, data))),
                false => None,
            })
    }

    /// Returns a parallel iterator on all top level nodes.
    /// # Examples
    /// ```
    /// use monstertruck_assembly::dag::*;
    /// use rayon::prelude::*;
    /// let mut dag = Dag::<usize, ()>::new();
    /// let a = dag.create_nodes(0..5);
    /// dag.create_edge(a[0], a[1], ());
    /// dag.create_edge(a[0], a[2], ());
    /// dag.create_edge(a[2], a[3], ());
    /// dag.create_edge(a[4], a[3], ());
    /// let tops = dag.par_top_nodes().collect::<Vec<_>>();
    /// assert_eq!(tops.len(), 2);
    /// assert_ne!(tops[0].entity(), tops[1].entity());
    /// ```
    #[inline]
    pub fn par_top_nodes(&self) -> ParTopNodes<'_, NE, EE>
    where
        NE: Sync,
        EE: Sync, {
        self.nodes
            .par_iter()
            .enumerate()
            .filter_map(move |(index, data)| match data.parents.is_empty() {
                true => Some(node((index, data))),
                false => None,
            })
    }

    /// Returns an iterator on all top level nodes.
    /// # Examples
    /// ```
    /// use monstertruck_assembly::dag::*;
    /// let mut dag = Dag::<usize, ()>::new();
    /// let a = dag.create_nodes(0..5);
    /// dag.create_edge(a[0], a[1], ());
    /// dag.create_edge(a[0], a[2], ());
    /// dag.create_edge(a[2], a[3], ());
    /// dag.create_edge(a[4], a[3], ());
    /// let tops = dag.top_nodes().collect::<Vec<_>>();
    /// assert_eq!(tops.len(), 2);
    /// assert_ne!(tops[0].entity(), tops[1].entity());
    /// ```
    #[inline]
    pub fn top_nodes_mut(&mut self) -> TopNodesMut<'_, NE, EE> {
        self.nodes
            .iter_mut()
            .enumerate()
            .filter_map(move |(index, data)| match data.parents.is_empty() {
                true => Some(node_mut((index, data))),
                false => None,
            })
    }

    /// Returns a parallel iterator on all top level nodes.
    /// # Examples
    /// ```
    /// use monstertruck_assembly::dag::*;
    /// use rayon::prelude::*;
    /// let mut dag = Dag::<usize, ()>::new();
    /// let a = dag.create_nodes(0..5);
    /// dag.create_edge(a[0], a[1], ());
    /// dag.create_edge(a[0], a[2], ());
    /// dag.create_edge(a[2], a[3], ());
    /// dag.create_edge(a[4], a[3], ());
    /// let tops = dag.par_top_nodes().collect::<Vec<_>>();
    /// assert_eq!(tops.len(), 2);
    /// assert_ne!(tops[0].entity(), tops[1].entity());
    /// ```
    #[inline]
    pub fn par_top_nodes_mut(&mut self) -> ParTopNodesMut<'_, NE, EE>
    where
        NE: Send,
        EE: Send, {
        self.nodes
            .par_iter_mut()
            .enumerate()
            .filter_map(move |(index, data)| match data.parents.is_empty() {
                true => Some(node_mut((index, data))),
                false => None,
            })
    }

    /// Add one node. Returns the index of the node.
    /// # Examples
    /// ```
    /// use monstertruck_assembly::dag::*;
    /// let mut dag = Dag::<usize, ()>::new();
    /// let index = dag.create_node(4);
    /// assert!(dag.node(index).is_some_and(|node| node.entity() == &4));
    /// ```
    pub fn create_node(&mut self, entity: NE) -> NodeIndex {
        let index = self.nodes.len();
        self.nodes.push(NodeData::new(entity));
        NodeIndex(index)
    }

    /// Add nodes from iterator. Returns the vector of the node indices.
    /// # Examples
    /// ```
    /// use monstertruck_assembly::dag::*;
    /// let mut dag = Dag::<usize, ()>::new();
    /// let a = dag.create_nodes(0..5);
    /// assert!(dag.node(a[3]).is_some_and(|node| node.entity() == &3));
    /// ```
    pub fn create_nodes<I>(&mut self, iter: I) -> Vec<NodeIndex>
    where I: IntoIterator<Item = NE> {
        let index = self.nodes.len();
        self.nodes.extend(iter.into_iter().map(NodeData::new));
        (index..self.nodes.len()).map(NodeIndex).collect()
    }

    /// Returns the node corresponding to `index`, or `None` if the index is
    /// out of range (for example, an index obtained from a larger graph).
    #[inline]
    pub fn node(&self, index: NodeIndex) -> Option<Node<'_, NE, EE>> {
        self.nodes.get(index.0).map(|data| Node { index, data })
    }

    /// Returns the mutable node corresponding to `index`, or `None` if the
    /// index is out of range.
    #[inline]
    pub fn node_mut(&mut self, index: NodeIndex) -> Option<NodeMut<'_, NE, EE>> {
        self.nodes
            .get_mut(index.0)
            .map(|data| NodeMut { index, data })
    }

    /// Add an edge. Returns `false` and do nothing if adding a edge causes a loop to occur.
    /// # Examples
    /// ```
    /// use monstertruck_assembly::dag::*;
    /// let mut dag = Dag::<(), usize>::new();
    /// let a = dag.create_nodes([(); 2]);
    /// assert!(dag.create_edge(a[0], a[1], 0));
    /// assert_eq!(dag.node(a[0]).map(|node| node.edges().len()), Some(1));
    ///
    /// // allow multiple edges between the same nodes
    /// assert!(dag.create_edge(a[0], a[1], 1));
    /// assert_eq!(dag.node(a[0]).map(|node| node.edges().len()), Some(2));
    ///
    /// // a loop has occurred!
    /// assert!(!dag.create_edge(a[1], a[0], 2));
    /// assert_eq!(dag.node(a[1]).map(|node| node.edges().len()), Some(0));
    /// ```
    pub fn create_edge(&mut self, from: NodeIndex, to: NodeIndex, entity: EE) -> bool {
        if self.node(from).is_none() || self.node(to).is_none() || self.has_path(to, from) {
            return false;
        }

        let Some(from_node) = self.node_mut(from) else {
            return false;
        };
        from_node.data.edges.push(EdgeData { to, entity });

        let Some(to_node) = self.node_mut(to) else {
            return false;
        };
        to_node.data.parents.insert(from);
        true
    }

    /// Remove an edge and returns its entity.
    /// # Examples
    /// ```
    /// use monstertruck_assembly::dag::*;
    /// let mut dag = Dag::<(), usize>::new();
    /// let a = dag.create_nodes([(); 2]);
    /// dag.create_edge(a[0], a[1], 0);
    /// dag.create_edge(a[0], a[1], 1);
    /// assert_eq!(dag.node(a[0]).map(|node| node.edges().len()), Some(2));
    ///
    /// let entity = dag.remove_edge(a[0], 1);
    /// assert_eq!(entity, Some(1));
    /// assert_eq!(dag.node(a[0]).map(|node| node.edges().len()), Some(1));
    /// ```
    pub fn remove_edge(&mut self, node: NodeIndex, index: usize) -> Option<EE> {
        let (to, entity, still_connected) = {
            let node_data = self.nodes.get_mut(node.0)?;
            if index >= node_data.edges.len() {
                return None;
            }
            let EdgeData { to, entity } = node_data.edges.remove(index);
            let still_connected = node_data.edges.iter().any(|edge_data| edge_data.to == to);
            (to, entity, still_connected)
        };
        if !still_connected && let Some(to_data) = self.nodes.get_mut(to.0) {
            to_data.parents.remove(&node);
        }
        Some(entity)
    }

    /// Returns an iterator on the all edges.
    /// # Examples
    /// ```
    /// use monstertruck_assembly::dag::*;
    /// let mut dag = Dag::<usize, ()>::new();
    /// let a = dag.create_nodes(0..5);
    /// dag.create_edge(a[0], a[1], ());
    /// dag.create_edge(a[0], a[2], ());
    /// dag.create_edge(a[2], a[3], ());
    /// dag.create_edge(a[4], a[3], ());
    /// let edges = dag.all_edges().collect::<Vec<_>>();
    /// assert_eq!(edges.len(), 4);
    /// ```
    #[inline]
    pub fn all_edges<'a>(
        &'a self,
    ) -> AllEdges<'a, NE, EE, impl FnMut(&'a EdgeData<EE>) -> Edge<'a, EE>> {
        self.all_nodes().flat_map(move |node| node.edges())
    }

    /// Returns a parallel iterator on the all edges.
    /// # Examples
    /// ```
    /// use monstertruck_assembly::dag::*;
    /// use rayon::prelude::*;
    /// let mut dag = Dag::<usize, ()>::new();
    /// let a = dag.create_nodes(0..5);
    /// dag.create_edge(a[0], a[1], ());
    /// dag.create_edge(a[0], a[2], ());
    /// dag.create_edge(a[2], a[3], ());
    /// dag.create_edge(a[4], a[3], ());
    /// let edges = dag.par_all_edges().collect::<Vec<_>>();
    /// assert_eq!(edges.len(), 4);
    /// ```
    #[inline]
    pub fn par_all_edges<'a>(
        &'a self,
    ) -> ParAllEdges<'a, NE, EE, impl Fn(&'a EdgeData<EE>) -> Edge<'a, EE>>
    where
        NE: Sync,
        EE: Sync, {
        self.par_all_nodes().flat_map(move |node| node.par_edges())
    }

    /// Returns an iterator for all paths with the top = `self`.
    /// # Examples
    /// ```
    /// use monstertruck_assembly::dag::*;
    /// let mut dag = Dag::<(), ()>::new();
    /// let a = dag.create_nodes([(); 5]);
    /// dag.create_edge(a[0], a[1], ());
    /// dag.create_edge(a[0], a[2], ());
    /// dag.create_edge(a[0], a[3], ());
    /// dag.create_edge(a[2], a[3], ());
    /// dag.create_edge(a[2], a[4], ());
    /// let slices = dag
    ///     .paths_iter(a[0])
    ///     .map(|path| {
    ///         path.nodes().iter().map(|node| node.index()).collect::<Vec<_>>()
    ///     })
    ///     .collect::<Vec<_>>();
    /// assert_eq!(slices.len(), 6);
    /// assert_eq!(slices[0], &[a[0]]);
    /// assert_eq!(slices[1], &[a[0], a[1]]);
    /// assert_eq!(slices[2], &[a[0], a[2]]);
    /// assert_eq!(slices[3], &[a[0], a[2], a[3]]);
    /// assert_eq!(slices[4], &[a[0], a[2], a[4]]);
    /// assert_eq!(slices[5], &[a[0], a[3]]);
    /// ```
    #[inline]
    pub fn paths_iter(&self, node: NodeIndex) -> impl Iterator<Item = Path<'_, NE, EE>> {
        PathsIter {
            dag: self,
            current: self.node(node).map(|node| PathIndices {
                nodes: vec![node],
                edges: Vec::new(),
            }),
            end_iter: false,
        }
    }

    /// Returns an iterator for all paths with the top = `self`.
    /// # Examples
    /// ```
    /// use monstertruck_assembly::dag::*;
    /// let mut dag = Dag::<(), ()>::new();
    /// let a = dag.create_nodes([(); 5]);
    /// dag.create_edge(a[0], a[1], ());
    /// dag.create_edge(a[0], a[2], ());
    /// dag.create_edge(a[0], a[3], ());
    /// dag.create_edge(a[2], a[3], ());
    /// dag.create_edge(a[2], a[4], ());
    /// let slices = dag
    ///     .paths_iter(a[0])
    ///     .map(|path| {
    ///         path.nodes().iter().map(|node| node.index()).collect::<Vec<_>>()
    ///     })
    ///     .collect::<Vec<_>>();
    /// assert_eq!(slices.len(), 6);
    /// assert_eq!(slices[0], &[a[0]]);
    /// assert_eq!(slices[1], &[a[0], a[1]]);
    /// assert_eq!(slices[2], &[a[0], a[2]]);
    /// assert_eq!(slices[3], &[a[0], a[2], a[3]]);
    /// assert_eq!(slices[4], &[a[0], a[2], a[4]]);
    /// assert_eq!(slices[5], &[a[0], a[3]]);
    /// ```
    #[inline]
    pub fn maximal_paths_iter(&self, node: NodeIndex) -> impl Iterator<Item = Path<'_, NE, EE>> {
        MaximalPathsIter {
            dag: self,
            current: self.first_maximal_path(node),
            end_iter: false,
        }
    }

    /// Returns true if a path exists from `from` to `to`.
    /// # Examples
    /// ```
    /// use monstertruck_assembly::dag::*;
    /// use rayon::prelude::*;
    /// let mut dag = Dag::<usize, ()>::new();
    /// let a = dag.create_nodes(0..5);
    /// dag.create_edge(a[0], a[1], ());
    /// dag.create_edge(a[0], a[2], ());
    /// dag.create_edge(a[2], a[3], ());
    /// dag.create_edge(a[4], a[3], ());
    /// assert!(!dag.has_path(a[0], a[4]));
    ///
    /// dag.create_edge(a[1], a[4], ());
    /// assert!(dag.has_path(a[0], a[4]));
    /// ```
    #[inline]
    pub fn has_path(&self, from: NodeIndex, to: NodeIndex) -> bool {
        if self.node(from).is_none() || self.node(to).is_none() {
            return false;
        }
        self.sub_has_path(from, to, &mut BTreeSet::default())
    }

    /// Returns true if a path exists from `from` to `to`.
    fn sub_has_path(&self, from: NodeIndex, to: NodeIndex, seen: &mut BTreeSet<NodeIndex>) -> bool {
        if from == to {
            return true;
        }
        let Some(to_node) = self.node(to) else {
            return false;
        };
        for &parent in to_node.parents() {
            if seen.contains(&parent) {
                continue;
            } else if self.sub_has_path(from, parent, seen) {
                return true;
            }
        }
        seen.insert(to);
        false
    }

    fn first_maximal_path(&self, from: NodeIndex) -> Option<PathIndices<'_, NE, EE>> {
        let mut cursor = self.node(from)?;
        let mut nodes = vec![cursor];
        let mut edges = vec![];
        while let Some(edge) = cursor.edge(0) {
            edges.push(0);
            cursor = self.node(edge.nodes().1)?;
            nodes.push(cursor);
        }
        Some(PathIndices { nodes, edges })
    }

    /// Creates a `Dag` based on node entities and adjacency lists.
    ///
    /// # Examples
    /// ```
    /// use monstertruck_assembly::dag::*;
    /// let adjacency = [
    ///     (0, 1, ()),
    ///     (0, 2, ()),
    ///     (1, 3, ()),
    ///     (1, 4, ()),
    ///     (4, 2, ()),
    /// ];
    ///
    /// let dag = Dag::try_from_adjacency(0..5, adjacency).unwrap();
    /// let nodes = dag.all_nodes().collect::<Vec<_>>();
    ///
    /// let parents = nodes[2].parents();
    /// assert_eq!(parents.len(), 2);
    /// assert!(parents.contains(&nodes[0].index()));
    /// assert!(parents.contains(&nodes[4].index()));
    /// ```
    ///
    /// # Failure
    /// Returns `None` if:
    /// - there exists an index in `adjacency` which is more than the length of `entities`, or,
    /// - the graph has a cycle.
    ///
    /// ```
    /// use monstertruck_assembly::dag::*;
    /// let adjacency = [
    ///     (0, 1, ()),
    ///     (0, 2, ()),
    ///     (1, 3, ()),
    ///     (1, 4, ()),
    ///     (4, 1, ()),
    ///     (4, 2, ()),
    /// ];
    /// assert!(Dag::try_from_adjacency(0..5, adjacency).is_none());
    /// ```
    pub fn try_from_adjacency(
        node_entities: impl IntoIterator<Item = NE>,
        adjacency: impl IntoIterator<Item = (usize, usize, EE)>,
    ) -> Option<Self> {
        let node_entities = node_entities.into_iter().collect::<Vec<_>>();
        let mut edges = node_entities
            .iter()
            .map(|_| Vec::<(usize, EE)>::new())
            .collect::<Vec<_>>();
        let mut parents = node_entities
            .iter()
            .map(|_| BTreeSet::<NodeIndex>::new())
            .collect::<Vec<_>>();
        for (i, j, edge_entity) in adjacency {
            if i >= node_entities.len() || j >= node_entities.len() {
                return None;
            }
            edges[i].push((j, edge_entity));
            parents[j].insert(NodeIndex(i));
        }
        if has_cycle(&parents) {
            return None;
        }
        let nodes = node_entities
            .into_iter()
            .zip(edges)
            .zip(parents)
            .map(|((node_entity, edges), parents)| NodeData {
                entity: node_entity,
                edges: edges
                    .into_iter()
                    .map(|(to, edge_entity)| EdgeData {
                        to: NodeIndex(to),
                        entity: edge_entity,
                    })
                    .collect(),
                parents,
            })
            .collect();
        Some(Self { nodes })
    }

    /// Create a new `Dag` by mapping node and edge entities to new values.
    /// # Examples
    /// ```
    /// use monstertruck_assembly::dag::*;
    /// let mut src_dag = Dag::<i32, i32>::new();
    /// let a = src_dag.create_nodes(0..3);
    /// src_dag.create_edge(a[0], a[1], 0);
    /// src_dag.create_edge(a[0], a[2], 1);
    /// src_dag.create_edge(a[1], a[2], 2);
    /// let dst_dag = src_dag.map(|n: &i32| 2f64.powi(*n), |n: &i32| 3f64.powi(*n));
    ///
    /// assert_eq!(dst_dag.len(), 3);
    /// let b = dst_dag.all_nodes().collect::<Vec<_>>();
    ///
    /// assert_eq!(b[0].index(), a[0]);
    /// assert_eq!(b[1].index(), a[1]);
    /// assert_eq!(b[2].index(), a[2]);
    ///
    /// assert_eq!(*b[0].entity(), 1.0);
    /// assert_eq!(*b[1].entity(), 2.0);
    /// assert_eq!(*b[2].entity(), 4.0);
    ///
    /// let edges0 = b[0].edges().collect::<Vec<_>>();
    /// assert_eq!(edges0.len(), 2);
    /// assert_eq!(edges0[0].nodes(), (b[0].index(), b[1].index()));
    /// assert_eq!(*edges0[0].entity(), 1.0);
    /// assert_eq!(edges0[1].nodes(), (b[0].index(), b[2].index()));
    /// assert_eq!(*edges0[1].entity(), 3.0);
    ///
    /// let edges1 = b[1].edges().collect::<Vec<_>>();
    /// assert_eq!(edges1.len(), 1);
    /// assert_eq!(edges1[0].nodes(), (b[1].index(), b[2].index()));
    /// assert_eq!(*edges1[0].entity(), 9.0);
    /// ```
    pub fn map<'a, NE2, EE2, NF, EF>(
        &'a self,
        mut node_map: NF,
        mut edge_map: EF,
    ) -> Dag<NE2, EE2>
    where
        NF: FnMut(&'a NE) -> NE2,
        EF: FnMut(&'a EE) -> EE2,
    {
        let mut edge_data_map = move |edge_data: &'a EdgeData<EE>| EdgeData {
            to: edge_data.to,
            entity: edge_map(&edge_data.entity),
        };
        let node_data_map = move |node_data: &'a NodeData<NE, EE>| NodeData {
            entity: node_map(&node_data.entity),
            edges: node_data.edges.iter().map(&mut edge_data_map).collect(),
            parents: node_data.parents.clone(),
        };

        Dag {
            nodes: self.nodes.iter().map(node_data_map).collect(),
        }
    }

    /// Create a new `Dag` by mapping node and edge entities to new values through parallel processing.
    /// # Examples
    /// ```
    /// use monstertruck_assembly::dag::*;
    /// let mut src_dag = Dag::<i32, i32>::new();
    /// let a = src_dag.create_nodes(0..3);
    /// src_dag.create_edge(a[0], a[1], 0);
    /// src_dag.create_edge(a[0], a[2], 1);
    /// src_dag.create_edge(a[1], a[2], 2);
    /// let dst_dag = src_dag.par_map(|n: &i32| 2f64.powi(*n), |n: &i32| 3f64.powi(*n));
    ///
    /// assert_eq!(dst_dag.len(), 3);
    /// let b = dst_dag.all_nodes().collect::<Vec<_>>();
    ///
    /// assert_eq!(b[0].index(), a[0]);
    /// assert_eq!(b[1].index(), a[1]);
    /// assert_eq!(b[2].index(), a[2]);
    ///
    /// assert_eq!(*b[0].entity(), 1.0);
    /// assert_eq!(*b[1].entity(), 2.0);
    /// assert_eq!(*b[2].entity(), 4.0);
    ///
    /// let edges0 = b[0].edges().collect::<Vec<_>>();
    /// assert_eq!(edges0.len(), 2);
    /// assert_eq!(edges0[0].nodes(), (b[0].index(), b[1].index()));
    /// assert_eq!(*edges0[0].entity(), 1.0);
    /// assert_eq!(edges0[1].nodes(), (b[0].index(), b[2].index()));
    /// assert_eq!(*edges0[1].entity(), 3.0);
    ///
    /// let edges1 = b[1].edges().collect::<Vec<_>>();
    /// assert_eq!(edges1.len(), 1);
    /// assert_eq!(edges1[0].nodes(), (b[1].index(), b[2].index()));
    /// assert_eq!(*edges1[0].entity(), 9.0);
    /// ```
    pub fn par_map<'a, NE2, EE2, NF, EF>(&'a self, node_map: NF, edge_map: EF) -> Dag<NE2, EE2>
    where
        NE: Sync,
        EE: Sync,
        NE2: Send,
        EE2: Send,
        NF: Fn(&'a NE) -> NE2 + Send + Sync,
        EF: Fn(&'a EE) -> EE2 + Send + Sync, {
        let edge_data_map = move |edge_data: &'a EdgeData<EE>| EdgeData {
            to: edge_data.to,
            entity: edge_map(&edge_data.entity),
        };
        let node_data_map = move |node_data: &'a NodeData<NE, EE>| NodeData {
            entity: node_map(&node_data.entity),
            edges: node_data.edges.iter().map(&edge_data_map).collect(),
            parents: node_data.parents.clone(),
        };

        Dag {
            nodes: self.nodes.par_iter().map(node_data_map).collect(),
        }
    }

    /// Create a new `Dag` by mapping node and edge entities to new values.
    /// # Examples
    /// ```
    /// use monstertruck_assembly::dag::*;
    /// let mut src_dag = Dag::<i32, i32>::new();
    /// let a = src_dag.create_nodes(0..3);
    /// src_dag.create_edge(a[0], a[1], 0);
    /// src_dag.create_edge(a[0], a[2], 1);
    /// src_dag.create_edge(a[1], a[2], 2);
    /// let dst_dag = src_dag.map_owned(|n: i32| 2f64.powi(n), |n: i32| 3f64.powi(n));
    ///
    /// assert_eq!(dst_dag.len(), 3);
    /// let b = dst_dag.all_nodes().collect::<Vec<_>>();
    ///
    /// assert_eq!(b[0].index(), a[0]);
    /// assert_eq!(b[1].index(), a[1]);
    /// assert_eq!(b[2].index(), a[2]);
    ///
    /// assert_eq!(*b[0].entity(), 1.0);
    /// assert_eq!(*b[1].entity(), 2.0);
    /// assert_eq!(*b[2].entity(), 4.0);
    ///
    /// let edges0 = b[0].edges().collect::<Vec<_>>();
    /// assert_eq!(edges0.len(), 2);
    /// assert_eq!(edges0[0].nodes(), (b[0].index(), b[1].index()));
    /// assert_eq!(*edges0[0].entity(), 1.0);
    /// assert_eq!(edges0[1].nodes(), (b[0].index(), b[2].index()));
    /// assert_eq!(*edges0[1].entity(), 3.0);
    ///
    /// let edges1 = b[1].edges().collect::<Vec<_>>();
    /// assert_eq!(edges1.len(), 1);
    /// assert_eq!(edges1[0].nodes(), (b[1].index(), b[2].index()));
    /// assert_eq!(*edges1[0].entity(), 9.0);
    /// ```
    pub fn map_owned<NE2, EE2, NF, EF>(self, mut node_map: NF, mut edge_map: EF) -> Dag<NE2, EE2>
    where
        NF: FnMut(NE) -> NE2,
        EF: FnMut(EE) -> EE2, {
        let mut edge_data_map = move |edge_data: EdgeData<EE>| EdgeData {
            to: edge_data.to,
            entity: edge_map(edge_data.entity),
        };
        let node_data_map = move |node_data: NodeData<NE, EE>| NodeData {
            entity: node_map(node_data.entity),
            edges: node_data
                .edges
                .into_iter()
                .map(&mut edge_data_map)
                .collect(),
            parents: node_data.parents.clone(),
        };

        Dag {
            nodes: self.nodes.into_iter().map(node_data_map).collect(),
        }
    }

    /// Create a new `Dag` by mapping node and edge entities to new values through parallel processing.
    /// # Examples
    /// ```
    /// use monstertruck_assembly::dag::*;
    /// let mut src_dag = Dag::<i32, i32>::new();
    /// let a = src_dag.create_nodes(0..3);
    /// src_dag.create_edge(a[0], a[1], 0);
    /// src_dag.create_edge(a[0], a[2], 1);
    /// src_dag.create_edge(a[1], a[2], 2);
    /// let dst_dag = src_dag.par_map_owned(|n: i32| 2f64.powi(n), |n: i32| 3f64.powi(n));
    ///
    /// assert_eq!(dst_dag.len(), 3);
    /// let b = dst_dag.all_nodes().collect::<Vec<_>>();
    ///
    /// assert_eq!(b[0].index(), a[0]);
    /// assert_eq!(b[1].index(), a[1]);
    /// assert_eq!(b[2].index(), a[2]);
    ///
    /// assert_eq!(*b[0].entity(), 1.0);
    /// assert_eq!(*b[1].entity(), 2.0);
    /// assert_eq!(*b[2].entity(), 4.0);
    ///
    /// let edges0 = b[0].edges().collect::<Vec<_>>();
    /// assert_eq!(edges0.len(), 2);
    /// assert_eq!(edges0[0].nodes(), (b[0].index(), b[1].index()));
    /// assert_eq!(*edges0[0].entity(), 1.0);
    /// assert_eq!(edges0[1].nodes(), (b[0].index(), b[2].index()));
    /// assert_eq!(*edges0[1].entity(), 3.0);
    ///
    /// let edges1 = b[1].edges().collect::<Vec<_>>();
    /// assert_eq!(edges1.len(), 1);
    /// assert_eq!(edges1[0].nodes(), (b[1].index(), b[2].index()));
    /// assert_eq!(*edges1[0].entity(), 9.0);
    /// ```
    pub fn par_map_owned<NE2, EE2, NF, EF>(self, node_map: NF, edge_map: EF) -> Dag<NE2, EE2>
    where
        NE: Send,
        EE: Send,
        NE2: Send,
        EE2: Send,
        NF: Fn(NE) -> NE2 + Send + Sync,
        EF: Fn(EE) -> EE2 + Send + Sync, {
        let edge_data_map = move |edge_data: EdgeData<EE>| EdgeData {
            to: edge_data.to,
            entity: edge_map(edge_data.entity),
        };
        let node_data_map = move |node_data: NodeData<NE, EE>| NodeData {
            entity: node_map(node_data.entity),
            edges: node_data.edges.into_iter().map(&edge_data_map).collect(),
            parents: node_data.parents.clone(),
        };

        Dag {
            nodes: self.nodes.into_par_iter().map(node_data_map).collect(),
        }
    }
}

fn has_cycle(adjacency: &[BTreeSet<NodeIndex>]) -> bool {
    let mut seen = vec![false; adjacency.len()];
    let mut finished = vec![false; adjacency.len()];

    while let Some(first) = finished.iter().position(|&x| !x) {
        if sub_has_cycle(first, adjacency, &mut seen, &mut finished) {
            return true;
        }
    }
    false
}

fn sub_has_cycle(
    cursor: usize,
    adjacency: &[BTreeSet<NodeIndex>],
    seen: &mut [bool],
    finished: &mut [bool],
) -> bool {
    seen[cursor] = true;
    for &NodeIndex(child) in &adjacency[cursor] {
        if finished[child] {
            continue;
        } else if seen[child] || sub_has_cycle(child, adjacency, seen, finished) {
            return true;
        }
    }
    finished[cursor] = true;
    false
}

impl<'a, NE, EE> Node<'a, NE, EE> {
    /// Returns node index
    #[inline]
    pub fn index(self) -> NodeIndex { self.index }
    /// Returns entity of node.
    #[inline]
    pub fn entity(self) -> &'a NE { &self.data.entity }
    /// Returns the `index`th edge.
    ///
    /// # Panics
    /// `index` must be smaller than the number of edges from `self`.
    #[inline]
    pub fn edge(self, index: usize) -> Option<Edge<'a, EE>> {
        Some(edge_closure(self.index)(self.data.edges.get(index)?))
    }
    /// Returns an iterator on edges from `self`.
    #[inline]
    pub fn edges(self) -> Edges<'a, EE, impl FnMut(&'a EdgeData<EE>) -> Edge<'a, EE>> {
        self.data.edges.iter().map(edge_closure(self.index))
    }
    /// Returns a parallel iterator on edges from `self`.
    #[inline]
    pub fn par_edges(self) -> ParEdges<'a, EE, impl Fn(&'a EdgeData<EE>) -> Edge<'a, EE>>
    where EE: Sync {
        self.data.edges.par_iter().map(edge_closure(self.index))
    }
    /// Returns parents from `self`.
    #[inline]
    pub fn parents(self) -> &'a BTreeSet<NodeIndex> { &self.data.parents }

    /// Returns `true` if `self` has no parents.
    /// # Examples
    /// ```
    /// use monstertruck_assembly::dag::*;
    /// let mut dag = Dag::<usize, usize>::new();
    /// let a = dag.create_nodes(0..2);
    /// dag.create_edge(a[0], a[1], 0);
    /// assert!(dag.node(a[0]).is_some_and(Node::is_top));
    /// assert!(dag.node(a[1]).is_some_and(|node| !node.is_top()));
    /// ```
    #[inline]
    pub fn is_top(self) -> bool { self.data.parents.is_empty() }

    /// Returns `true` if `self` has no child.
    /// # Examples
    /// ```
    /// use monstertruck_assembly::dag::*;
    /// let mut dag = Dag::<usize, usize>::new();
    /// let a = dag.create_nodes(0..2);
    /// dag.create_edge(a[0], a[1], 0);
    /// assert!(dag.node(a[1]).is_some_and(Node::is_terminal));
    /// assert!(dag.node(a[0]).is_some_and(|node| !node.is_terminal()));
    /// ```
    #[inline]
    pub fn is_terminal(self) -> bool { self.data.edges.is_empty() }
}

impl<'a, NE, EE> NodeMut<'a, NE, EE> {
    /// Returns node index
    #[inline]
    pub fn index(&self) -> NodeIndex { self.index }
    /// Returns entity of node.
    #[inline]
    pub fn entity(&mut self) -> &mut NE { &mut self.data.entity }
    /// Returns the `index`th edge.
    pub fn edge(&'a mut self, index: usize) -> Option<EdgeMut<'a, EE>> {
        Some(edge_mut_closure(self.index)(
            self.data.edges.get_mut(index)?,
        ))
    }
    /// Returns an iterator on edges from `self`.
    #[inline]
    pub fn edges(&'a self) -> Edges<'a, EE, impl FnMut(&'a EdgeData<EE>) -> Edge<'a, EE>> {
        self.data.edges.iter().map(edge_closure(self.index))
    }
    /// Returns an iterator on edges from `self`.
    #[inline]
    pub fn edges_mut(
        &'a mut self,
    ) -> EdgesMut<'a, EE, impl FnMut(&'a mut EdgeData<EE>) -> EdgeMut<'a, EE>> {
        self.data.edges.iter_mut().map(edge_mut_closure(self.index))
    }
    /// Returns a parallel iterator on edges from `self`.
    #[inline]
    pub fn par_edges(&'a self) -> ParEdges<'a, EE, impl Fn(&'a EdgeData<EE>) -> Edge<'a, EE>>
    where EE: Sync {
        self.data.edges.par_iter().map(edge_closure(self.index))
    }
    /// Returns a parallel iterator on edges from `self`.
    #[inline]
    pub fn par_edges_mut(
        &'a mut self,
    ) -> ParEdgesMut<'a, EE, impl Fn(&'a mut EdgeData<EE>) -> EdgeMut<'a, EE>>
    where EE: Send {
        self.data
            .edges
            .par_iter_mut()
            .map(edge_mut_closure(self.index))
    }
    /// Returns parents from `self`.
    #[inline]
    pub fn parents(&'a self) -> &'a BTreeSet<NodeIndex> { &self.data.parents }

    /// Returns `true` if `self` has no parents.
    #[inline]
    pub fn is_top(&self) -> bool { self.data.parents.is_empty() }

    /// Returns `true` if `self` has no child.
    #[inline]
    pub fn is_terminal(&self) -> bool { self.data.edges.is_empty() }
}

impl<'a, EE> Edge<'a, EE> {
    /// Returns end nodes
    /// # Examples
    /// ```
    /// use monstertruck_assembly::dag::*;
    /// let mut dag = Dag::<usize, usize>::new();
    /// let a = dag.create_nodes(0..2);
    /// dag.create_edge(a[0], a[1], 0);
    /// let Some(edge) = dag.node(a[0]).and_then(|node| node.edge(0)) else { return };
    /// assert_eq!(edge.nodes(), (a[0], a[1]));
    /// ```
    #[inline]
    pub fn nodes(self) -> (NodeIndex, NodeIndex) { self.nodes }
    /// Returns entity
    /// # Examples
    /// ```
    /// use monstertruck_assembly::dag::*;
    /// let mut dag = Dag::<usize, usize>::new();
    /// let a = dag.create_nodes(0..2);
    /// dag.create_edge(a[0], a[1], 3);
    /// let Some(edge) = dag.node(a[0]).and_then(|node| node.edge(0)) else { return };
    /// assert_eq!(edge.entity(), &3);
    /// ```
    #[inline]
    pub fn entity(self) -> &'a EE { self.entity }
}

impl<'a, EE> EdgeMut<'a, EE> {
    /// Returns end nodes
    /// # Examples
    /// ```
    /// use monstertruck_assembly::dag::*;
    /// let mut dag = Dag::<usize, usize>::new();
    /// let a = dag.create_nodes(0..2);
    /// dag.create_edge(a[0], a[1], 0);
    /// let Some(mut node) = dag.node_mut(a[0]) else { return };
    /// let Some(edge) = node.edge(0) else { return };
    /// assert_eq!(edge.nodes(), (a[0], a[1]));
    /// ```
    #[inline]
    pub fn nodes(&self) -> (NodeIndex, NodeIndex) { self.nodes }
    /// Returns entity
    /// # Examples
    /// ```
    /// use monstertruck_assembly::dag::*;
    /// let mut dag = Dag::<usize, usize>::new();
    /// let a = dag.create_nodes(0..2);
    /// dag.create_edge(a[0], a[1], 3);
    /// {
    ///     let Some(mut node) = dag.node_mut(a[0]) else { return };
    ///     let Some(mut edge) = node.edge(0) else { return };
    ///     let entity = edge.entity();
    ///     assert_eq!(entity, &3);
    ///     *entity = 5;
    /// }
    /// assert!(dag
    ///     .node(a[0])
    ///     .and_then(|node| node.edge(0))
    ///     .is_some_and(|edge| edge.entity() == &5));
    /// ```
    #[inline]
    pub fn entity(&'a mut self) -> &'a mut EE { self.entity }
}

impl<'a, NE, EE> Path<'a, NE, EE> {
    /// Returns nodes of path.
    #[inline]
    pub fn nodes(&self) -> &[Node<'a, NE, EE>] { &self.nodes }
    /// Returns an iterator of node indices.
    #[inline]
    pub fn node_indices(&self) -> impl Iterator<Item = NodeIndex> {
        self.nodes.iter().map(|node| node.index)
    }
    /// Returns nodes of path.
    #[inline]
    pub fn edges(&self) -> &[Edge<'a, EE>] { &self.edges }
    /// Returns the terminal node
    #[inline]
    pub fn terminal_node(&'a self) -> Option<Node<'a, NE, EE>> { self.nodes.last().copied() }
}

#[derive(Debug)]
struct PathIndices<'a, NE, EE> {
    nodes: Vec<Node<'a, NE, EE>>,
    edges: Vec<usize>,
}

impl<'a, NE, EE> Clone for PathIndices<'a, NE, EE> {
    fn clone(&self) -> Self {
        Self {
            nodes: self.nodes.clone(),
            edges: self.edges.clone(),
        }
    }
}

impl<'a, NE, EE> PathIndices<'a, NE, EE> {
    fn into_path(self) -> Option<Path<'a, NE, EE>> {
        let nodes = self.nodes;
        let edges = self
            .edges
            .into_iter()
            .zip(&nodes)
            .map(|(index, node)| node.edge(index))
            .collect::<Option<Vec<_>>>()?;
        Some(Path { nodes, edges })
    }
}

#[derive(Clone, Debug)]
struct PathsIter<'a, NE, EE> {
    dag: &'a Dag<NE, EE>,
    current: Option<PathIndices<'a, NE, EE>>,
    end_iter: bool,
}

impl<'a, NE, EE> Iterator for PathsIter<'a, NE, EE> {
    type Item = Path<'a, NE, EE>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.end_iter {
            return None;
        }
        let current = self.current.as_mut()?;
        let res = current.clone();
        let leaf = *current.nodes.last()?;

        if leaf.is_terminal() {
            let some_idx = current
                .edges
                .iter()
                .enumerate()
                .rev()
                .find_map(|(idx, edge_idx)| {
                    let node = current.nodes.get(idx)?;
                    (*edge_idx + 1 < node.edges().len()).then_some(idx)
                });
            if let Some(idx) = some_idx {
                current.nodes.truncate(idx + 2);
                current.edges.truncate(idx + 1);
                let edge_idx = current.edges.get_mut(idx)?;
                *edge_idx += 1;
                let node = *current.nodes.get(idx)?;
                let (_, to) = node.edge(*edge_idx)?.nodes();
                let next = self.dag.node(to)?;
                *current.nodes.get_mut(idx + 1)? = next;
            } else {
                self.end_iter = true;
            }
        } else {
            let (_, new_leaf_index) = leaf.edge(0)?.nodes();
            current.nodes.push(self.dag.node(new_leaf_index)?);
            current.edges.push(0);
        }

        let result = res.into_path();
        if result.is_none() {
            self.end_iter = true;
        }
        result
    }
}

#[derive(Clone, Debug)]
struct MaximalPathsIter<'a, NE, EE> {
    dag: &'a Dag<NE, EE>,
    current: Option<PathIndices<'a, NE, EE>>,
    end_iter: bool,
}

impl<'a, NE, EE> Iterator for MaximalPathsIter<'a, NE, EE> {
    type Item = Path<'a, NE, EE>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.end_iter {
            return None;
        }
        let current = self.current.as_mut()?;
        let res = current.clone();

        let some_idx = current
            .edges
            .iter()
            .enumerate()
            .rev()
            .find_map(|(idx, edge_idx)| {
                let node = current.nodes.get(idx)?;
                (*edge_idx + 1 < node.edges().len()).then_some(idx)
            });
        let Some(idx) = some_idx else {
            self.end_iter = true;
            return res.into_path();
        };

        current.nodes.truncate(idx + 1);
        current.edges.truncate(idx + 1);
        let edge_idx = current.edges.get_mut(idx)?;
        *edge_idx += 1;

        let current_leaf = *current.nodes.get(idx)?;
        let (_, next) = current_leaf.edge(*edge_idx)?.nodes();
        let path = self.dag.first_maximal_path(next)?;
        current.nodes.extend(path.nodes);
        current.edges.extend(path.edges.into_iter().map(|_| 0));

        let result = res.into_path();
        if result.is_none() {
            self.end_iter = true;
        }
        result
    }
}

type NodeIndices = std::iter::Map<std::ops::Range<usize>, fn(usize) -> NodeIndex>;
type ParNodeIndices = rayon::iter::Map<rayon::range::Iter<usize>, fn(usize) -> NodeIndex>;

type AllNodes<'a, NE, EE> = std::iter::Map<
    std::iter::Enumerate<std::slice::Iter<'a, NodeData<NE, EE>>>,
    fn((usize, &NodeData<NE, EE>)) -> Node<'_, NE, EE>,
>;

type AllNodesMut<'a, NE, EE> = std::iter::Map<
    std::iter::Enumerate<std::slice::IterMut<'a, NodeData<NE, EE>>>,
    fn((usize, &mut NodeData<NE, EE>)) -> NodeMut<'_, NE, EE>,
>;

type ParAllNodes<'a, NE, EE> = rayon::iter::Map<
    rayon::iter::Enumerate<rayon::slice::Iter<'a, NodeData<NE, EE>>>,
    fn((usize, &NodeData<NE, EE>)) -> Node<'_, NE, EE>,
>;

type ParAllNodesMut<'a, NE, EE> = rayon::iter::Map<
    rayon::iter::Enumerate<rayon::slice::IterMut<'a, NodeData<NE, EE>>>,
    fn((usize, &mut NodeData<NE, EE>)) -> NodeMut<'_, NE, EE>,
>;

type TopNodes<'a, NE, EE> = std::iter::FilterMap<
    std::iter::Enumerate<std::slice::Iter<'a, NodeData<NE, EE>>>,
    fn((usize, &NodeData<NE, EE>)) -> Option<Node<'_, NE, EE>>,
>;

type TopNodesMut<'a, NE, EE> = std::iter::FilterMap<
    std::iter::Enumerate<std::slice::IterMut<'a, NodeData<NE, EE>>>,
    fn((usize, &mut NodeData<NE, EE>)) -> Option<NodeMut<'_, NE, EE>>,
>;

type ParTopNodes<'a, NE, EE> = rayon::iter::FilterMap<
    rayon::iter::Enumerate<rayon::slice::Iter<'a, NodeData<NE, EE>>>,
    fn((usize, &NodeData<NE, EE>)) -> Option<Node<'_, NE, EE>>,
>;

type ParTopNodesMut<'a, NE, EE> = rayon::iter::FilterMap<
    rayon::iter::Enumerate<rayon::slice::IterMut<'a, NodeData<NE, EE>>>,
    fn((usize, &mut NodeData<NE, EE>)) -> Option<NodeMut<'_, NE, EE>>,
>;

type Edges<'a, EE, F> = std::iter::Map<std::slice::Iter<'a, EdgeData<EE>>, F>;
type ParEdges<'a, EE, F> = rayon::iter::Map<rayon::slice::Iter<'a, EdgeData<EE>>, F>;
type EdgesMut<'a, EE, F> = std::iter::Map<std::slice::IterMut<'a, EdgeData<EE>>, F>;
type ParEdgesMut<'a, EE, F> = rayon::iter::Map<rayon::slice::IterMut<'a, EdgeData<EE>>, F>;

type AllEdges<'a, NE, EE, F> = std::iter::FlatMap<
    AllNodes<'a, NE, EE>,
    Edges<'a, EE, F>,
    fn(Node<'a, NE, EE>) -> Edges<'a, EE, F>,
>;

type ParAllEdges<'a, NE, EE, F> =
    rayon::iter::FlatMap<ParAllNodes<'a, NE, EE>, fn(Node<'a, NE, EE>) -> ParEdges<'a, EE, F>>;
