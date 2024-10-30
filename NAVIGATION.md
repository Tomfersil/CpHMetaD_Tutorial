```mermaid
flowchart LR
A[Metadynamics Intro] ==> B[CpH-Metadynamics intro]
B ==> C[CpH-MetaD source]
B ==> D[Conformational Analysis]
E ==> E[Protonation Analysis]
D ==> F[Solution]

click A "ref1" "A brief introduction to metadynamics is useful before executing this tutorial"
click B "README.md" "A theoretical introduction to the CpH-Metadynamics method"
click B "ref3" "Take a look into the source page and material for the CpH-Metadynamics"
click D "Conformation.md" "Perform the conformation analysis of CpH-Metadynamics."
click E "Protonation.md" "Perform the protonation analysis of CpH-Metadynamics"
click F "notebooks/solution.ipynb" "A python notebook containing a full set of solutions for the exercise"
```
