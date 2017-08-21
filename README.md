# nonOrthogonalityCheckSalome
A python class to be included in Salome meshing scripts for detection of non-orthogonal and skew faces.

Scheme:
MeshQualityCheck is checking your mesh for non-orthogonal and skew faces. The algorithm loops through all mesh elements and their respective neighbors to find internal faces. This is a long-lasting but necessary task as internal faces are not defined in Salome meshes. For the non-orthogonality test each internal face the angle between the line connecting the two face-neighboring elements' center of mass ('COMs line') and the face's normal direction is calculated and compared to the threshold. For the skewness test the distance between the COMs line and the face's center of mass is calculated, normalized with the length of the COMS line and compared to the skewness threshold. If a face is failing either of these tests, its neighboring elements are grouped and added to the list of failed elements.

Usage:
Include MeshQualityCheck to your meshing script. Make an instance of that class for every computed mesh part that you want to find the non-orthogonal and skew faces of. To define the the MeshQualityCheck instance you provide the mesh part and the thresholds for non-orthogonality (default value = 65) and skewness (default value = 0.5) as arguments. To deactivate the search for non-orthogonal or skew faces you set its respective threshold to 0. The search is started by calling checkMesh on the MeshQualityCheck instance. The pairs of mesh elements that share a face which failed the check are grouped and can be easily visualized in Salome's object viewer.

Example Case:
Find the attached case mashCheckTestCase.py for an example of MeshQualityCheck's application: A wedge geometry is created and meshed including a viscous layer on the wall surface. The mesh is computed, exported as a unv file (for comparison with OpenFOAM only) and finally checked for non-orthogonal and skew faces.

Motive:
Salome is a great tool for geometry and mesh generation for OpenFOAM simulations. But because of OpenFOAM's sensitivity of the mesh quality, and especially the non-orthogonality and skewness of its faces, I want to find bad faces or elements already in the meshing stage, prior to the mesh conversion. 


