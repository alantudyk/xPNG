#include "../until_fork/until_fork.h"

MAIN_VOID {
    
    ssize_t
     Y_min = 1000,  Y_max = -1000,
    Co_min = 1000, Co_max = -1000,
    Cg_min = 1000, Cg_max = -1000;
    
    fiN(R, 256)
        fiN(G, 256)
            fiN(B, 256) {
                            
                ssize_t Y, Co, Cg, r, g, b;
                
                Cg = G - (Y = B + ((Co = R - B) >> 1)), Y += Cg >> 1;
                
                if (Y < Y_min) Y_min = Y;
                if (Y > Y_max) Y_max = Y;
                if (Co < Co_min) Co_min = Co;
                if (Co > Co_max) Co_max = Co;
                if (Cg < Cg_min) Cg_min = Cg;
                if (Cg > Cg_max) Cg_max = Cg;
                
                g = Cg + (b = Y - (Cg >> 1)), r = (b -= Co >> 1) + Co;
                
                if (R != r || G != g || B != b) {
                    pf("\nFailed: R = %4ld,  G = %4ld,  B = %4ld\n",   R,  G,  B);
                    pf(  "        r = %4ld,  g = %4ld,  b = %4ld\n",   r,  g,  b);
                    pf(  "        Y = %4ld, Co = %4ld, Cg = %4ld\n\n", Y, Co, Cg);
                    ret 1;
                }
                
            }
    
    pf("\n"
    
    " Y_min = %4ld,  Y_max = %4ld,\n"
    "Co_min = %4ld, Co_max = %4ld,\n"
    "Cg_min = %4ld, Cg_max = %4ld.\n"
    
    , Y_min, Y_max, Co_min, Co_max, Cg_min, Cg_max);
    
    puts("\nOK.\n");
    
    ret 0;
}
