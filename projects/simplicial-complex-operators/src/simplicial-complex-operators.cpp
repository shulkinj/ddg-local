// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"
#include "geometrycentral/numerical/linear_algebra_utilities.h"
#include <Eigen/Sparse>
#include <cassert>


using namespace geometrycentral;
using namespace geometrycentral::surface;

/*
 * Assign a unique index to each vertex, edge, and face of a mesh.
 * All elements are 0-indexed.
 *
 * Input: None. Access geometry via the member variable <geometry>, and pointer to the mesh via <mesh>.
 * Returns: None.
 */
void SimplicialComplexOperators::assignElementIndices() {

    // Needed to access geometry->vertexIndices, etc. as cached quantities.
    // Not needed if you're just using v->getIndex(), etc.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    // You can set the index field of a vertex via geometry->vertexIndices[v], where v is a Vertex object (or an
    // integer). Similarly you can do edges and faces via geometry->edgeIndices, geometry->faceIndices, like so:
    size_t idx = 0;
    for (Vertex v : mesh->vertices()) {
        idx = geometry->vertexIndices[v];
    }

    for (Edge e : mesh->edges()) {
        idx = geometry->edgeIndices[e];
    }

    for (Face f : mesh->faces()) {
        idx = geometry->faceIndices[f];
    }

    // You can more easily get the indices of mesh elements using the function getIndex(), albeit less efficiently and
    // technically less safe (although you don't need to worry about it), like so:
    //
    //      v.getIndex()
    //
    // where v can be a Vertex, Edge, Face, Halfedge, etc. For example:

    for (Vertex v : mesh->vertices()) {
        idx = v.getIndex(); // == geometry->vertexIndices[v])
    }

    // Geometry Central already sets the indices for us, though, so this function is just here for demonstration.
    // You don't have to do anything :)
}

/*
 * Construct the unsigned vertex-edge adjacency matrix A0.
 *
 * Input:
 * Returns: The sparse vertex-edge adjacency matrix which gets stored in the global variable A0.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {

    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();
    // Note: You can build an Eigen sparse matrix from triplets, then return it as a Geometry Central SparseMatrix.
    // See <https://eigen.tuxfamily.org/dox/group__TutorialSparse.html> for documentation.
    
    // initialize eigensparse matrix
    SparseMatrix<size_t> mat(mesh->nEdges(),mesh->nVertices());

    for (Edge e : mesh->edges()){
        // we compute by edge since we know there will only be 2 per row
        // enter sparse matrix values one by one:
        size_t i1 = geometry->vertexIndices[e.firstVertex()];
        size_t i2 = geometry->vertexIndices[e.secondVertex()];
        size_t e_i = geometry->edgeIndices[e];
 
        mat.insert(e_i,i1)=1;
        mat.insert(e_i,i2)=1;

    }

    return mat; 
}

/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input:
 * Returns: The sparse face-edge adjacency matrix which gets stored in the global variable A1.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {
    
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    // initialize eigensparse matrix
    SparseMatrix<size_t> mat(mesh->nFaces(),mesh->nEdges());
    
    for (Halfedge he : mesh->halfedges()){

        // only count interior edges
        if(he.isInterior()){
        // enter sparse matrix values one by one:
        size_t e = geometry->edgeIndices[he.edge()];
        
        size_t f = geometry->faceIndices[he.face()];
        
        mat.insert(f,e)=1;
        }
        
    }
    return mat; 
}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {
    
    Vector<size_t> vertices= Vector<size_t>::Zero(mesh->nVertices(),1);

    std::set<size_t>::iterator it;
    
    for (it = subset.vertices.begin(); it!= subset.vertices.end(); it++){
        vertices[*it]=1;
    }
    
    return vertices;
}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {

    Vector<size_t> edges= Vector<size_t>::Zero(mesh->nEdges(),1);

    std::set<size_t>::iterator it;
    
    for (it = subset.edges.begin(); it!= subset.edges.end(); it++){
         edges[*it]=1;
    }
    
    return edges;
}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {

    Vector<size_t> faces= Vector<size_t>::Zero(mesh->nFaces(),1);

    std::set<size_t>::iterator it;
    
    for (it = subset.faces.begin(); it!= subset.faces.end(); it++){
         faces[*it]=1;
    }
    
    return faces;
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 */
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {


    MeshSubset star = subset;
    
    Vector<size_t> E0;
    Vector<size_t> F0;
    Vector<size_t> F1;

    std::set<size_t> set_E0;
    std::set<size_t> set_F0;
    std::set<size_t> set_F1;


    Vector<size_t> sub_verts = buildVertexVector(subset);
    Vector<size_t> sub_edges = buildEdgeVector(subset);
    Vector<size_t> sub_faces = buildFaceVector(subset);
   
    // Compute all connected vertices in the vectors
    E0 = A0 * sub_verts;
    F0 = A1 * A0 * sub_verts;
    F1 = A1 * sub_edges;


    // fill edge sets
    for(int i =0; i<E0.rows(); ++i){
        if(E0[i] != 0){
            set_E0.insert(i);
        }
    } 
    // fill face sets
    for(int i =0; i<F0.rows(); ++i){
        if(F0[i] != 0){
            set_F0.insert(i);
        }
        if(F1[i] != 0){
            set_F1.insert(i);
        }
    }
    
    star.addEdges(set_E0);
    star.addFaces(set_F0);
    star.addFaces(set_F1);

    return star;
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {
    
    // final returned subset
    MeshSubset closure = subset;
    
    // vectors for computation
    Vector<size_t> V0;
    Vector<size_t> E1;
    Vector<size_t> V1;

    // sets for storing indices
    std::set<size_t> set_V0;
    std::set<size_t> set_E1;
    std::set<size_t> set_V1;

    // Only need the subset edges and faces to compute the closure
    Vector<size_t> sub_edges = buildEdgeVector(subset);
    Vector<size_t> sub_faces = buildFaceVector(subset);
    
    // Add vertices that touch edges
    V0 = A0.transpose()*sub_edges;
    for(int i =0; i<V0.rows(); ++i){
        if(V0[i] != 0){
            set_V0.insert(i);
        }
    }
    closure.addVertices(set_V0);


    // Add edges that touch faces
    E1 = A1.transpose()*sub_faces;
    for(int i =0; i<E1.rows(); i++){
        if(E1[i] != 0){
            set_E1.insert(i);
        }
    }
    closure.addEdges(set_E1);

    // Add edges that touch faces
    V1 =A0.transpose()*A1.transpose()*sub_faces;
    for(int i =0; i<V1.rows(); i++){
        if(V1[i] != 0){
            set_V1.insert(i);
        }
    }
    closure.addVertices(set_V1);


    return closure;
    }

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {

    // define Cl(St(subset)) and St(Cl(subset))
    MeshSubset ClSt = closure(star(subset));
    MeshSubset StCl = star(closure(subset));

    // to remove StCl from ClSt we pick out the set of faces, edges and vertices of StCl
    std::set<size_t> v = StCl.vertices;
    std::set<size_t> e = StCl.edges;
    std::set<size_t> f = StCl.faces;
    
    // initialize link as Cl(St(subset))
    MeshSubset link = ClSt;
    
    link.deleteVertices(v);
    link.deleteEdges(e);
    link.deleteFaces(f);

    return link;
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {
    
    MeshSubset Cl;
    Cl = closure(subset);
    
    // a potentially faster method that only compares sizes rather than each elt
    size_t v0 = subset.vertices.size();
    size_t e0 = subset.edges.size();
    size_t f0 = subset.faces.size();
    size_t v1 = Cl.vertices.size();
    size_t e1 = Cl.edges.size();
    size_t f1 = Cl.faces.size();

    return (v0==v1 && e0==e1 && f0==f1);
    //return closure==subset;

}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */
int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {

    if(isComplex(subset)){
        Vector<size_t> sub_verts = buildVertexVector(subset);
        Vector<size_t> sub_edges = buildEdgeVector(subset);
        Vector<size_t> sub_faces = buildFaceVector(subset);

        
        if(subset.faces.size()!=0){
            // make sure every edge and vertex is part of a face
            Vector<size_t> new_edges = A1.transpose() * sub_faces;
            for(int i =0; i<new_edges.rows(); i++){
                if(new_edges[i] != 0){
                    new_edges[i]=1;
                }
            }
            Vector<size_t> new_verts = A0.transpose() * A1.transpose() * sub_faces;
            for(int i =0; i<new_verts.rows(); i++){
                if(new_verts[i] != 0){
                    new_verts[i]=1; 
                }
            }

            if(new_edges == sub_edges && new_verts == sub_verts){
                return 2;
            }
            else{
                return -1;   
            }
        }

        else if(subset.edges.size()!=0){
            // make sure every vertex is part of an edge
            Vector<size_t> new_verts = A0.transpose() * sub_edges;
            for(int i =0; i<new_verts.rows(); i++){
                if(new_verts[i] != 0){
                    new_verts[i]=1;
                }
            }

            if(new_verts == sub_verts){
                return 1;
            }
            else{
                return -1;
            }
        }

        // if only vertices we are done
        else{
            return 0;
        }
    }
    return -1;
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {
    
    MeshSubset boundary;

    // find simplicial complex degree
    int deg = isPureComplex(subset);
    
    // if not pure complex, return empty placeholder
    // (or error?)
    if(deg == -1){
        return boundary;
    }
    
    // if subset is points, return empty placeholder
    else if(deg == 0){
        return boundary;    
    }
    
    // 1-simplicial complex, return points connected to one edge
    else if(deg ==1){
        // only vertices that touch 1 edge
        // first find vertex connection vector
        Vector<size_t> edges = buildEdgeVector(subset);
        Vector<size_t> vertices = A0.transpose()*edges;
        std::set<size_t> int_verts_set;
        
        
        // fill bd_edges_set and turn edges into boundary edges 
        for(int i =0; i<vertices.rows(); i++){
            if(vertices[i] == 1){
                int_verts_set.insert(i);
            }
        }
        boundary.addVertices(int_verts_set);
    }
    
    // since boundaries are simplicial complexes of degree one lower:
    // we can only find the boundary edges, then take the closure at the end
    else if(deg ==2){ 
        // first find edges connected to one face
        // then take the closure
        Vector<size_t> faces = buildFaceVector(subset);
        Vector<size_t> edges = A1.transpose()*faces;
        std::set<size_t> int_edges_set;
        MeshSubset bdry_only_edges;
        
        
        for(int i =0; i<edges.rows(); i++){
            if(edges[i] == 1){
                int_edges_set.insert(i);
            }
        }
        bdry_only_edges.addEdges(int_edges_set);

        boundary = closure(bdry_only_edges);
        
    }
        
    return boundary; 
    }
