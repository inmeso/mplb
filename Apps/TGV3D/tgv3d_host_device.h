
#ifndef TGV3D_HOST_DEVICE_H
#define TGV3D_HOST_DEVICE_H
#ifndef OPS_FUN_PREFIX
#define OPS_FUN_PREFIX
#endif
#define Central2nd(p1, m1, dx) (((p1) - (m1)) / (2 * (dx)));

#define Central4th(p1, p2, m1, m2, dx) \
    ((-(p2) + 8 * (p1)-8 * (m1) + (m2)) / (12 * (dx)));

#define Central6th(p1, p2, p3, m1, m2, m3, dx) \
    (((p3)-9 * (p2) + 45 * (p1)-45 * (m1) + 9 * (m2) - (m3)) / (60 * (dx)));

#define Central8th(p1, p2, p3, p4, m1, m2, m3, m4, dx)                  \
    ((-(p4) / 280) + (4 * (p3) / 105) + (-(p2) / 5) + (4 * (p1) / 5) +  \
     ((m4) / 280) + (-4 * (m3) / 105) + ((m2) / 5) + (-4 * (m1) / 5)) / \
        (dx);

#endif //  TGV3D_HOST_DEVICE_H