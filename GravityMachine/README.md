
## Note on computing the object's displacement relative to the fluid.

### Prequisites: 
- Stage velocity relative to lab.
- Object's velocity with respect to the lab.
- Fluid velocity with respect to the lab.

Note that in gravity machine v2.0 and higher, the optical system is fixed relative to lab.
Therefore displacements relative to the optical system are displacements relative to the lab.
For earlier Gravity machine versions the optical system is fixed in Z but moving in X and Y. 

Stage velocity relative to lab (recorded in track data: Z)
Fluid velocity relative to lab (from PIV)
Object's velocity relative to lab (from object motion in the image)

We want to compute the object's velocity relative to the fluid.

Vector operations to calculate V_objFluid which is what we want.

Note: Vz is the object velocity relative to the stage V_objStage
V_objLab: is measured from the displacement of the object centroid in the image. 
V_objStage = V_objLab - V_stageLab
Therefore, 
 V_stageLab = V_objLab - VobjStage   ---- (1)

The measured fluid velocity using PIV is:
 V_measured = V_stageLab + V_fluidStage
 Therefore, 
 V_fluidStage = V_measured - V_stageLab ---- (2)
 We can substitute for V_stageLab from (1) to get V_fluidStage
 
Now, 
	 V_objFluid = V_objStage - V_fluidStage 
	