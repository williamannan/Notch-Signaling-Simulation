# Notch-Signaling-Simulation
The code simulates the Notch-Delta interaction model with filopodia dynamics included

The outline for the code is as follows:

 "TwoDGeom.m" - create a 2D tissue packed with hexagonal cells. It stores the center coordinates as well as the coordinates at the vertices of the cells.

 "Filop_vectors.m" - This function file computes the directional vectors along the filopodia in the tissue, the coordinates at the base and tips of the filopodia, and the apical and filopodia neighbors for each cell in the tissue

"ParDeltaIN.m"- This function computes the effective Delta expressed by the neighboring cells. It performs the computation in a parallel manner

"FilopAngles.m" - This function computes the angle distribution of the filopodia when oriented towards a fixed point in the tissue. Here, we selected five fixed locations, namely, northeast (NE), northwest (NW), southwest (SW), southeast (SE) and center (C). To orient all filopodia to a fixed point, you need to pass "Ang" as an input parameter to the "GenTrendFilopDyn.m" function. You also have to comment the "Ang" defined in this function. In the output file, you need to uncomment this line "[Ang] = FilopAngles(Nc, r, 'C')". 

"FilopLent.m" - This function solves the filopodia length dynamic model

"GenTrendFilopDyn.m"- This is the main function file that solves the Notch-Delta signaling model

"OutPut.m" -  This is where all the model parameters are specified and also the various function files are called to run.

 


