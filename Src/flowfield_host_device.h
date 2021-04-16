#ifndef FLOWFIELD_HOST_DEVICE_H
#define FLOWFIELD_HOST_DEVICE_H
#ifndef OPS_FUN_PREFIX
#define OPS_FUN_PREFIX
#endif
#define NEWVERTEX
#ifdef NEWVERTEX
// Notes on periodic boundary condition
// There could be two types of implementation, one is same to finite difference
// scheme, the other one is same to the molecular dynamics.
// We provide both methods.
// Finite difference way:
// Using the kernel function  KerCutCellPeriodic3D when treating the boundary
// condition. A halo transfer for g_f is required before treating boundary
// condition. In this case, the grid is treated as if a wall node
// Molecular dynamics:
// In this way, nothing is needed when treating the boundary condition. However,
// a halo transfer is needed for g_fStage before the stream step following the
// current implementation stream-collision. In this case, the grid is treated
// as if a fluid node
enum class VertexType {
    // vtfluid is the general type of node
    // All specific fluid types should be started as vtft
    Fluid = 10,
    Inlet = 11,
    OutLet = 12,
    MDPeriodic = 13,
    FDPeriodic = 15,
    Symmetry = 14,
    Wall = 1000,
    // Bulk node but with immersed solid nodes
    ImmersedSolid = -1,
    // Immersed solid but surface
    ImmersedBoundary = -2,
    //surface and/or edge between two blocks
    VirtualBoundary = -3
};  // vt
#endif
#ifndef NEWVERTEX
enum VertexType {
    // vtfluid is the general type of node
    // All specific fluid types should be started as vtft
    VertexType::Fluid = 10,
    Vertex_Boundary = 1000,
    // Vertex_Boundary is the general type of node
    // All specific boundary types should be started as vtbt
    Vertex_KineticDiffuseWall = 1001,
    Vertex_KineticSpelluarWall = 1002,
    Vertex_SlipWall = 1003,
    Vertex_VelocityInlet = 1004,
    Vertex_VelocityOutlet = 1005,
    Vertex_ExtrapolPressure1ST = 1006,
    Vertex_ExtrapolPressure2ND = 1007,
    Vertex_Periodic = 1008,
    Vertex_Uniform = 1009,
    Vertex_BounceBackWall = 1010,
    Vertex_FreeFlux = 1011,
    Vertex_ZouHeVelocity = 1012,
    Vertex_EQMDiffuseRefl = 1014,
    Vertex_NoslipEQN = 1020,
    // Vertex_BoundaryCorner is the general corner type of node
    // All specific corner types should be started as vtbtco
    // for example, a corner of bounceback wall and inlet may be named as
    // vtbtbbinlet
    // It is unlikely to be used.
    Vertex_BoundaryCorner = 2000,
    // Bulk node but with immersed solid nodes
    VertexType::ImmersedSolid = -1
};  // vt
#endif
// Use this type to describe the geometry property of a node in terms of surface
// for example, it is i=0, j=0, or i=imax,j=jmax
// Rule i(x)=1 j(y)=2 k(z)=3, start(min)=0 end(max)=1
// So the surface i=0 is 10
// in the 2D case, we suppose that there are only i,j planes.
//
enum VertexGeometryType {
    // boundary node which usually a envelop of the region
    // boundary surface
    // VertexGeometry= VG
    // a normal fluid node, no need to specify geometry property
    VG_Fluid = 0,
    // a normal immersed solid node, no need to specify geometry property
    VG_ImmersedSolid = -1,
    // surface property, the normal direction
    // the normal direction pointing to the Positive (P) X(I) direction
    VG_IP = 10,
    // the normal direction pointing to the Minor (M) X(I) direction
    VG_IM = 11,
    VG_JP = 20,
    VG_JM = 21,
    VG_KP = 30,  // 3D
    VG_KM = 31,  // 3D
    // for inner corners
    // boundary cross line for 3D and corner point for 2D
    // only meaningful for 3D lines
    VG_IPKP_I = 10300,
    VG_IPKM_I = 10310,
    VG_IMKP_I = 11300,
    VG_IMKM_I = 11310,
    VG_JPKP_I = 20300,
    VG_JPKM_I = 20310,
    VG_JMKP_I = 21300,
    VG_JMKM_I = 21310,
    // for both 3D lines and 2D points
    VG_IPJP_I = 10200,
    VG_IPJM_I = 10210,
    VG_IMJP_I = 11200,
    VG_IMJM_I = 11210,
    // boundary corner point for 3D
    VG_IPJPKP_I = 1020300,
    VG_IPJPKM_I = 1020310,
    VG_IPJMKP_I = 1021300,
    VG_IPJMKM_I = 1021310,
    VG_IMJPKP_I = 1120300,
    VG_IMJPKM_I = 1120310,
    VG_IMJMKP_I = 1121300,
    VG_IMJMKM_I = 1121310,
    // for outer corners
    // only meaningful for 3D lines
    VG_IPKP_O = 10301,
    VG_IPKM_O = 10311,
    VG_IMKP_O = 11301,
    VG_IMKM_O = 11311,
    VG_JPKP_O = 20301,
    VG_JPKM_O = 20311,
    VG_JMKP_O = 21301,
    VG_JMKM_O = 21311,
    // for both 3D lines and 2D points
    VG_IPJP_O = 10201,
    VG_IPJM_O = 10211,
    VG_IMJP_O = 11201,
    VG_IMJM_O = 11211,
    // boundary corner point for 3D
    VG_IPJPKP_O = 1020301,
    VG_IPJPKM_O = 1020311,
    VG_IPJMKP_O = 1021301,
    VG_IPJMKM_O = 1021311,
    VG_IMJPKP_O = 1120301,
    VG_IMJPKM_O = 1120311,
    VG_IMJMKP_O = 1121301,
    VG_IMJMKM_O = 1121311,

};  // vg
#ifdef OPS_3D
static inline OPS_FUN_PREFIX int SpaceDim(){return 3;};
#endif
#ifdef OPS_2D
static inline OPS_FUN_PREFIX int SpaceDim(){return 2;};
#endif

#endif // FLOWFIELD_HOST_DEVICE_H