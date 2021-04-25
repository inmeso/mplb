#ifndef BOUNDARY_HOST_DEVICE_H
#define BOUNDARY_HOST_DEVICE_H
#ifndef OPS_FUN_PREFIX
#define OPS_FUN_PREFIX
#endif

enum class BoundarySurface {
    Left = 0,
    Right = 1,
    Top = 2,
    Bottom = 3,
#ifdef OPS_3D
    Front = 4,
    Back = 5,
    LeftBack = 50,
    LeftFront = 40,
    RightBack = 51,
    RightFront = 41,
    TopBack = 52,
    TopFront = 42,
    BottomBack = 53,
    BottomFront = 43,
#endif
    LeftTop = 20,
    LeftBottom = 30,
    RightTop = 21,
    RightBottom = 31

};

/*!
 * Discrete velocity type at a solid wall boundary
 */
enum BndryDvType {
    BndryDv_Incoming = 1,
    BndryDv_Outgoing = 2,
    BndryDv_Parallel = 3,
    BndryDv_Invalid = -1
};
/*!
 * @brief Determining discrete velocity type a solid wall boundary
 * @param vg Geometry property, e.g., corner type
 * @param discreteVelocity components of a discrete velocity
 * @return BndryDvType discrete velocity type
 */

static inline OPS_FUN_PREFIX  BndryDvType FindBdyDvType(const VertexGeometryType vg,
                          const Real* discreteVelocity) {
    Real cx = discreteVelocity[0];
    Real cy = discreteVelocity[1];
    BndryDvType res{BndryDv_Invalid};
    switch (vg) {
        case VG_IP: {
            if (cx < 0) {
                res = BndryDv_Incoming;
            }
            if (cx == 0) {
                res = BndryDv_Parallel;
            }
            if (cx > 0) {
                res = BndryDv_Outgoing;
            }
        } break;
        case VG_IM: {
            if (cx > 0) {
                res = BndryDv_Incoming;
            }
            if (cx == 0) {
                res = BndryDv_Parallel;
            }
            if (cx < 0) {
                res = BndryDv_Outgoing;
            }
        } break;
        case VG_JP: {
            if (cy < 0) {
                res = BndryDv_Incoming;
            }
            if (cy == 0) {
                res = BndryDv_Parallel;
            }
            if (cy > 0) {
                res = BndryDv_Outgoing;
            }

        } break;
        case VG_JM: {
            if (cy > 0) {
                res = BndryDv_Incoming;
            }
            if (cy == 0) {
                res = BndryDv_Parallel;
            }
            if (cy < 0) {
                res = BndryDv_Outgoing;
            }
        } break;
        case VG_IPJP_I: {
            if ((cx <= 0 && cy < 0) || (cx < 0 && cy == 0)) {
                res = BndryDv_Incoming;
            } else if ((cx >= 0 && cy > 0) || (cx > 0 && cy == 0)) {
                res = BndryDv_Outgoing;
            } else {
                res = BndryDv_Parallel;
            }
        } break;
        case VG_IPJM_I: {
            if ((cx <= 0 && cy > 0) || (cx < 0 && cy == 0)) {
                res = BndryDv_Incoming;
            } else if ((cx >= 0 && cy < 0) || (cx > 0 && cy == 0)) {
                res = BndryDv_Outgoing;
            } else {
                res = BndryDv_Parallel;
            }
        } break;
        case VG_IMJP_I: {
            if ((cx >= 0 && cy < 0) || (cx > 0 && cy == 0)) {
                res = BndryDv_Incoming;
            } else if ((cx <= 0 && cy > 0) || (cx < 0 && cy == 0)) {
                res = BndryDv_Outgoing;
            } else {
                res = BndryDv_Parallel;
            }
        } break;
        case VG_IMJM_I: {
            if ((cx >= 0 && cy > 0) || (cx > 0 && cy == 0)) {
                res = BndryDv_Incoming;
            } else if ((cx <= 0 && cy < 0) || (cx < 0 && cy == 0)) {
                res = BndryDv_Outgoing;
            } else {
                res = BndryDv_Parallel;
            }

        } break;
        case VG_IPJP_O: {
            if ((cx <= 0 && cy < 0) || (cx < 0 && cy == 0)) {
                res = BndryDv_Incoming;
            } else if ((cx >= 0 && cy > 0) || (cx > 0 && cy == 0)) {
                res = BndryDv_Outgoing;
            } else {
                res = BndryDv_Parallel;
            }
        } break;
        case VG_IPJM_O: {
            if ((cx <= 0 && cy > 0) || (cx < 0 && cy == 0)) {
                res = BndryDv_Incoming;
            } else if ((cx >= 0 && cy < 0) || (cx > 0 && cy == 0)) {
                res = BndryDv_Outgoing;
            } else {
                res = BndryDv_Parallel;
            }
        } break;
        case VG_IMJP_O: {
            if ((cx >= 0 && cy < 0) || (cx > 0 && cy == 0)) {
                res = BndryDv_Incoming;
            } else if ((cx <= 0 && cy > 0) || (cx < 0 && cy == 0)) {
                res = BndryDv_Outgoing;
            } else {
                res = BndryDv_Parallel;
            }
        } break;
        case VG_IMJM_O: {
            if ((cx >= 0 && cy > 0) || (cx > 0 && cy == 0)) {
                res = BndryDv_Incoming;
            } else if ((cx <= 0 && cy < 0) || (cx < 0 && cy == 0)) {
                res = BndryDv_Outgoing;
            } else {
                res = BndryDv_Parallel;
            }

        } break;
        default:
            break;
    }
#ifdef debug
    if (BndryDv_Invalid == res) {
        ops_printf("%s\n",
                   "Warning: Some wrong in determining discrete velocity type "
                   "at 2D boundary ");
    }
#endif
    return res;
}

/*!
 * @brief Determining discrete velocity type a solid wall boundary: 3D
 * @param vg Geometry property, e.g., corner type
 * @param discreteVelocity components of a discrete velocity
 * @return BndryDvType discrete velocity type
 */
static inline OPS_FUN_PREFIX BndryDvType FindBdyDvType3D(const VertexGeometryType vg,
                            const Real* discreteVelocity) {
    Real cx = discreteVelocity[0];
    Real cy = discreteVelocity[1];
    Real cz = discreteVelocity[2];
    BndryDvType res{BndryDv_Invalid};
    switch (vg) {
        case VG_IP: {
            if (cx < 0) {
                res = BndryDv_Incoming;
            }
            if (cx == 0) {
                res = BndryDv_Parallel;
            }
            if (cx > 0) {
                res = BndryDv_Outgoing;
            }
        } break;
        case VG_IM: {
            if (cx > 0) {
                res = BndryDv_Incoming;
            }
            if (cx == 0) {
                res = BndryDv_Parallel;
            }
            if (cx < 0) {
                res = BndryDv_Outgoing;
            }
        } break;
        case VG_JP: {
            if (cy < 0) {
                res = BndryDv_Incoming;
            }
            if (cy == 0) {
                res = BndryDv_Parallel;
            }
            if (cy > 0) {
                res = BndryDv_Outgoing;
            }

        } break;
        case VG_JM: {
            if (cy > 0) {
                res = BndryDv_Incoming;
            }
            if (cy == 0) {
                res = BndryDv_Parallel;
            }
            if (cy < 0) {
                res = BndryDv_Outgoing;
            }
        } break;
        case VG_KP: {
            if (cz < 0) {
                res = BndryDv_Incoming;
            }
            if (cz == 0) {
                res = BndryDv_Parallel;
            }
            if (cz > 0) {
                res = BndryDv_Outgoing;
            }

        } break;
        case VG_KM: {
            if (cz > 0) {
                res = BndryDv_Incoming;
            }
            if (cz == 0) {
                res = BndryDv_Parallel;
            }
            if (cz < 0) {
                res = BndryDv_Outgoing;
            }
        } break;
        case VG_IPJP_I: {
            if ((cx <= 0 && cy < 0) || (cx < 0 && cy == 0)) {
                res = BndryDv_Incoming;
            } else if ((cx >= 0 && cy > 0) || (cx > 0 && cy == 0)) {
                res = BndryDv_Outgoing;
            } else {
                res = BndryDv_Parallel;
            }
        } break;
        case VG_IPJM_I: {
            if ((cx <= 0 && cy > 0) || (cx < 0 && cy == 0)) {
                res = BndryDv_Incoming;
            } else if ((cx >= 0 && cy < 0) || (cx > 0 && cy == 0)) {
                res = BndryDv_Outgoing;
            } else {
                res = BndryDv_Parallel;
            }
        } break;
        case VG_IMJP_I: {
            if ((cx >= 0 && cy < 0) || (cx > 0 && cy == 0)) {
                res = BndryDv_Incoming;
            } else if ((cx <= 0 && cy > 0) || (cx < 0 && cy == 0)) {
                res = BndryDv_Outgoing;
            } else {
                res = BndryDv_Parallel;
            }
        } break;
        case VG_IMJM_I: {
            if ((cx >= 0 && cy > 0) || (cx > 0 && cy == 0)) {
                res = BndryDv_Incoming;
            } else if ((cx <= 0 && cy < 0) || (cx < 0 && cy == 0)) {
                res = BndryDv_Outgoing;
            } else {
                res = BndryDv_Parallel;
            }

        } break;
        case VG_IPKP_I: {
            if ((cx <= 0 && cz < 0) || (cx < 0 && cz == 0)) {
                res = BndryDv_Incoming;
            } else if ((cx >= 0 && cz > 0) || (cx > 0 && cz == 0)) {
                res = BndryDv_Outgoing;
            } else {
                res = BndryDv_Parallel;
            }
        } break;
        case VG_IPKM_I: {
            if ((cx <= 0 && cz > 0) || (cx < 0 && cz == 0)) {
                res = BndryDv_Incoming;
            } else if ((cx >= 0 && cz < 0) || (cx > 0 && cz == 0)) {
                res = BndryDv_Outgoing;
            } else {
                res = BndryDv_Parallel;
            }
        } break;
        case VG_IMKP_I: {
            if ((cx >= 0 && cz < 0) || (cx > 0 && cz == 0)) {
                res = BndryDv_Incoming;
            } else if ((cx <= 0 && cz > 0) || (cx < 0 && cz == 0)) {
                res = BndryDv_Outgoing;
            } else {
                res = BndryDv_Parallel;
            }
        } break;
        case VG_IMKM_I: {
            if ((cx >= 0 && cz > 0) || (cx > 0 && cz == 0)) {
                res = BndryDv_Incoming;
            } else if ((cx <= 0 && cz < 0) || (cx < 0 && cz == 0)) {
                res = BndryDv_Outgoing;
            } else {
                res = BndryDv_Parallel;
            }

        } break;
        case VG_JPKP_I: {
            if ((cy <= 0 && cz < 0) || (cy < 0 && cz == 0)) {
                res = BndryDv_Incoming;
            } else if ((cy >= 0 && cz > 0) || (cy > 0 && cz == 0)) {
                res = BndryDv_Outgoing;
            } else {
                res = BndryDv_Parallel;
            }
        } break;
        case VG_JPKM_I: {
            if ((cy <= 0 && cz > 0) || (cy < 0 && cz == 0)) {
                res = BndryDv_Incoming;
            } else if ((cy >= 0 && cz < 0) || (cy > 0 && cz == 0)) {
                res = BndryDv_Outgoing;
            } else {
                res = BndryDv_Parallel;
            }
        } break;
        case VG_JMKP_I: {
            if ((cy >= 0 && cz < 0) || (cy > 0 && cz == 0)) {
                res = BndryDv_Incoming;
            } else if ((cy <= 0 && cz > 0) || (cy < 0 && cz == 0)) {
                res = BndryDv_Outgoing;
            } else {
                res = BndryDv_Parallel;
            }
        } break;
        case VG_JMKM_I: {
            if ((cy >= 0 && cz > 0) || (cy > 0 && cz == 0)) {
                res = BndryDv_Incoming;
            } else if ((cy <= 0 && cz < 0) || (cy < 0 && cz == 0)) {
                res = BndryDv_Outgoing;
            } else {
                res = BndryDv_Parallel;
            }

        } break;
        case VG_IPJP_O: {
            if ((cx <= 0 && cy < 0) || (cx < 0 && cy == 0)) {
                res = BndryDv_Incoming;
            } else if ((cx >= 0 && cy > 0) || (cx > 0 && cy == 0)) {
                res = BndryDv_Outgoing;
            } else {
                res = BndryDv_Parallel;
            }
        } break;
        case VG_IPJM_O: {
            if ((cx <= 0 && cy > 0) || (cx < 0 && cy == 0)) {
                res = BndryDv_Incoming;
            } else if ((cx >= 0 && cy < 0) || (cx > 0 && cy == 0)) {
                res = BndryDv_Outgoing;
            } else {
                res = BndryDv_Parallel;
            }
        } break;
        case VG_IMJP_O: {
            if ((cx >= 0 && cy < 0) || (cx > 0 && cy == 0)) {
                res = BndryDv_Incoming;
            } else if ((cx <= 0 && cy > 0) || (cx < 0 && cy == 0)) {
                res = BndryDv_Outgoing;
            } else {
                res = BndryDv_Parallel;
            }
        } break;
        case VG_IMJM_O: {
            if ((cx >= 0 && cy > 0) || (cx > 0 && cy == 0)) {
                res = BndryDv_Incoming;
            } else if ((cx <= 0 && cy < 0) || (cx < 0 && cy == 0)) {
                res = BndryDv_Outgoing;
            } else {
                res = BndryDv_Parallel;
            }

        } break;
        case VG_IPKP_O: {
            if ((cx <= 0 && cz < 0) || (cx < 0 && cz == 0)) {
                res = BndryDv_Incoming;
            } else if ((cx >= 0 && cz > 0) || (cx > 0 && cz == 0)) {
                res = BndryDv_Outgoing;
            } else {
                res = BndryDv_Parallel;
            }
        } break;
        case VG_IPKM_O: {
            if ((cx <= 0 && cz > 0) || (cx < 0 && cz == 0)) {
                res = BndryDv_Incoming;
            } else if ((cx >= 0 && cz < 0) || (cx > 0 && cz == 0)) {
                res = BndryDv_Outgoing;
            } else {
                res = BndryDv_Parallel;
            }
        } break;
        case VG_IMKP_O: {
            if ((cx >= 0 && cz < 0) || (cx > 0 && cz == 0)) {
                res = BndryDv_Incoming;
            } else if ((cx <= 0 && cz > 0) || (cx < 0 && cz == 0)) {
                res = BndryDv_Outgoing;
            } else {
                res = BndryDv_Parallel;
            }
        } break;
        case VG_IMKM_O: {
            if ((cx >= 0 && cz > 0) || (cx > 0 && cz == 0)) {
                res = BndryDv_Incoming;
            } else if ((cx <= 0 && cz < 0) || (cx < 0 && cz == 0)) {
                res = BndryDv_Outgoing;
            } else {
                res = BndryDv_Parallel;
            }
        } break;
        case VG_JPKP_O: {
            if ((cx <= 0 && cz < 0) || (cx < 0 && cz == 0)) {
                res = BndryDv_Incoming;
            } else if ((cx >= 0 && cz > 0) || (cx > 0 && cz == 0)) {
                res = BndryDv_Outgoing;
            } else {
                res = BndryDv_Parallel;
            }
        } break;
        case VG_JPKM_O: {
            if ((cy <= 0 && cz > 0) || (cy < 0 && cz == 0)) {
                res = BndryDv_Incoming;
            } else if ((cy >= 0 && cz < 0) || (cy > 0 && cz == 0)) {
                res = BndryDv_Outgoing;
            } else {
                res = BndryDv_Parallel;
            }
        } break;
        case VG_JMKP_O: {
            if ((cy >= 0 && cz < 0) || (cy > 0 && cz == 0)) {
                res = BndryDv_Incoming;
            } else if ((cy <= 0 && cz > 0) || (cy < 0 && cz == 0)) {
                res = BndryDv_Outgoing;
            } else {
                res = BndryDv_Parallel;
            }
        } break;
        case VG_JMKM_O: {
            if ((cy >= 0 && cz > 0) || (cy > 0 && cz == 0)) {
                res = BndryDv_Incoming;
            } else if ((cy <= 0 && cz < 0) || (cy < 0 && cz == 0)) {
                res = BndryDv_Outgoing;
            } else {
                res = BndryDv_Parallel;
            }
        } break;
        case VG_IPJPKP_I: {
            if (cx == 0 && cy == 0 && cz == 0) {
                res = BndryDv_Parallel;
            } else {
                if (cx <= 0 && cy <= 0 && cz <= 0) {
                    res = BndryDv_Outgoing;
                } else if (cx >= 0 && cy >= 0 && cz >= 0) {
                    res = BndryDv_Incoming;
                } else {
                    res = BndryDv_Parallel;
                }
            }
        } break;
        case VG_IPJPKM_I: {
            if (cx == 0 && cy == 0 && cz == 0) {
                res = BndryDv_Parallel;
            } else {
                if (cx <= 0 && cy <= 0 && cz >= 0) {
                    res = BndryDv_Outgoing;
                } else if (cx >= 0 && cy >= 0 && cz <= 0) {
                    res = BndryDv_Incoming;
                } else {
                    res = BndryDv_Parallel;
                }
            }
        } break;
        case VG_IPJMKP_I: {
            if (cx == 0 && cy == 0 && cz == 0) {
                res = BndryDv_Parallel;
            } else {
                if (cx <= 0 && cy >= 0 && cz <= 0) {
                    res = BndryDv_Outgoing;
                } else if (cx >= 0 && cy <= 0 && cz >= 0) {
                    res = BndryDv_Incoming;
                } else {
                    res = BndryDv_Parallel;
                }
            }
        } break;
        case VG_IPJMKM_I: {
            if (cx == 0 && cy == 0 && cz == 0) {
                res = BndryDv_Parallel;
            } else {
                if (cx <= 0 && cy >= 0 && cz >= 0) {
                    res = BndryDv_Outgoing;
                } else if (cx >= 0 && cy <= 0 && cz <= 0) {
                    res = BndryDv_Incoming;
                } else {
                    res = BndryDv_Parallel;
                }
            }
        } break;
        case VG_IMJPKP_I: {
            if (cx == 0 && cy == 0 && cz == 0) {
                res = BndryDv_Parallel;
            } else {
                if (cx >= 0 && cy <= 0 && cz <= 0) {
                    res = BndryDv_Outgoing;
                } else if (cx <= 0 && cy >= 0 && cz >= 0) {
                    res = BndryDv_Incoming;
                } else {
                    res = BndryDv_Parallel;
                }
            }
        } break;
        case VG_IMJPKM_I: {
            if (cx == 0 && cy == 0 && cz == 0) {
                res = BndryDv_Parallel;
            } else {
                if (cx >= 0 && cy <= 0 && cz >= 0) {
                    res = BndryDv_Outgoing;
                } else if (cx <= 0 && cy >= 0 && cz <= 0) {
                    res = BndryDv_Incoming;
                } else {
                    res = BndryDv_Parallel;
                }
            }
        } break;
        case VG_IMJMKP_I: {
            if (cx == 0 && cy == 0 && cz == 0) {
                res = BndryDv_Parallel;
            } else {
                if (cx >= 0 && cy >= 0 && cz <= 0) {
                    res = BndryDv_Outgoing;
                } else if (cx <= 0 && cy <= 0 && cz >= 0) {
                    res = BndryDv_Incoming;
                } else {
                    res = BndryDv_Parallel;
                }
            }
        } break;
        case VG_IMJMKM_I: {
            if (cx == 0 && cy == 0 && cz == 0) {
                res = BndryDv_Parallel;
            } else {
                if (cx >= 0 && cy >= 0 && cz >= 0) {
                    res = BndryDv_Outgoing;
                } else if (cx <= 0 && cy <= 0 && cz <= 0) {
                    res = BndryDv_Incoming;
                } else {
                    res = BndryDv_Parallel;
                }
            }
        } break;
        case VG_IPJPKP_O: {
            if (cx == 0 && cy == 0 && cz == 0) {
                res = BndryDv_Parallel;
            } else {
                if (cx <= 0 && cy <= 0 && cz <= 0) {
                    res = BndryDv_Outgoing;
                } else if (cx >= 0 && cy >= 0 && cz >= 0) {
                    res = BndryDv_Incoming;
                } else {
                    res = BndryDv_Parallel;
                }
            }
        } break;
        case VG_IPJPKM_O: {
            if (cx == 0 && cy == 0 && cz == 0) {
                res = BndryDv_Parallel;
            } else {
                if (cx <= 0 && cy <= 0 && cz >= 0) {
                    res = BndryDv_Outgoing;
                } else if (cx >= 0 && cy >= 0 && cz <= 0) {
                    res = BndryDv_Incoming;
                } else {
                    res = BndryDv_Parallel;
                }
            }
        } break;
        case VG_IPJMKP_O: {
            if (cx == 0 && cy == 0 && cz == 0) {
                res = BndryDv_Parallel;
            } else {
                if (cx <= 0 && cy >= 0 && cz <= 0) {
                    res = BndryDv_Outgoing;
                } else if (cx >= 0 && cy <= 0 && cz >= 0) {
                    res = BndryDv_Incoming;
                } else {
                    res = BndryDv_Parallel;
                }
            }
        } break;
        case VG_IPJMKM_O: {
            if (cx == 0 && cy == 0 && cz == 0) {
                res = BndryDv_Parallel;
            } else {
                if (cx <= 0 && cy >= 0 && cz >= 0) {
                    res = BndryDv_Outgoing;
                } else if (cx >= 0 && cy <= 0 && cz <= 0) {
                    res = BndryDv_Incoming;
                } else {
                    res = BndryDv_Parallel;
                }
            }
        } break;
        case VG_IMJPKP_O: {
            if (cx == 0 && cy == 0 && cz == 0) {
                res = BndryDv_Parallel;
            } else {
                if (cx >= 0 && cy <= 0 && cz <= 0) {
                    res = BndryDv_Outgoing;
                } else if (cx <= 0 && cy >= 0 && cz >= 0) {
                    res = BndryDv_Incoming;
                } else {
                    res = BndryDv_Parallel;
                }
            }
        } break;
        case VG_IMJPKM_O: {
            if (cx == 0 && cy == 0 && cz == 0) {
                res = BndryDv_Parallel;
            } else {
                if (cx >= 0 && cy <= 0 && cz >= 0) {
                    res = BndryDv_Outgoing;
                } else if (cx <= 0 && cy >= 0 && cz <= 0) {
                    res = BndryDv_Incoming;
                } else {
                    res = BndryDv_Parallel;
                }
            }
        } break;
        case VG_IMJMKP_O: {
            if (cx == 0 && cy == 0 && cz == 0) {
                res = BndryDv_Parallel;
            } else {
                if (cx >= 0 && cy >= 0 && cz <= 0) {
                    res = BndryDv_Outgoing;
                } else if (cx <= 0 && cy <= 0 && cz >= 0) {
                    res = BndryDv_Incoming;
                } else {
                    res = BndryDv_Parallel;
                }
            }
        } break;
        case VG_IMJMKM_O: {
            if (cx == 0 && cy == 0 && cz == 0) {
                res = BndryDv_Parallel;
            } else {
                if (cx >= 0 && cy >= 0 && cz >= 0) {
                    res = BndryDv_Outgoing;
                } else if (cx <= 0 && cy <= 0 && cz <= 0) {
                    res = BndryDv_Incoming;
                } else {
                    res = BndryDv_Parallel;
                }
            }
        } break;

        default:
            break;
    }
#ifdef debug
    if (BndryDv_Invalid == res) {
        ops_printf("%s\n",
                   "Warning: Something wrong in determining discrete velocity "
                   "type at 3D boundary");
    }
#endif
    return res;
}
#endif //  BOUNDARY_HOST_DEVICE_H