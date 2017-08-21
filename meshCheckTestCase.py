# meshCheckTestCase.py
# Copyright (C) 2017  Fynn Aschmoneit

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


import salome
salome.salome_init()
import GEOM, SMESH

from salome.geom import geomBuilder
geompy = geomBuilder.New(salome.myStudy)

from salome.smesh import smeshBuilder
smesh =  smeshBuilder.New(salome.myStudy)

from meshQualityCheck import MeshQualityCheck

p1 = geompy.MakeVertex(0, 0, 0)
p2 = geompy.MakeVertex(5, 0, 2)
p3 = geompy.MakeVertex(5, 0, 0)

line1 = geompy.MakeLineTwoPnt( p1, p2 )
line2 = geompy.MakeLineTwoPnt( p1, p3 )
line3 = geompy.MakeLineTwoPnt( p2, p3 )

face = geompy.MakeFace( geompy.MakeWire( [line1, line2, line3] ), True )
vol = geompy.MakePrismDXDYDZ(face, 0, 6, 0, theScaleFactor = -1.0)
partition = geompy.MakePartition( [vol] )
faceList = geompy.SubShapeAllSorted(partition, geompy.ShapeType["FACE"])
for i in range(0,len(faceList)):
	geompy.addToStudy(faceList[i], "face_no_" + str(i) )

mesh = smesh.Mesh( partition )

algo1D = mesh.Segment()
algo2D = mesh.Triangle()
algo3D = mesh.Tetrahedron()

EBL = 1		# egde base length

algo1D.LocalLength(EBL)
algo3D.ViscousLayers(EBL*0.2, 1, 1.6, [ faceList[2], faceList[3] ] )

grDict = {}

grDict[	"inlet" ] 	= geompy.CreateGroup( partition, geompy.ShapeType["FACE"] )
grDict[	"outlet" ] 	= geompy.CreateGroup( partition, geompy.ShapeType["FACE"] )
grDict[	"wall" ] 	= geompy.CreateGroup( partition, geompy.ShapeType["FACE"] )

geompy.AddObject( grDict[	"inlet" ] , geompy.GetSubShapeID(partition, faceList[2]) )
geompy.AddObject( grDict[	"outlet" ] , geompy.GetSubShapeID(partition, faceList[3]) )
geompy.AddObject( grDict[	"wall" ] , geompy.GetSubShapeID(partition, faceList[0] ) )
geompy.AddObject( grDict[	"wall" ] , geompy.GetSubShapeID(partition, faceList[1] ) )
geompy.AddObject( grDict[	"wall" ] , geompy.GetSubShapeID(partition, faceList[4] ) )

mesh.Compute()

for k, v in grDict.iteritems():
		mesh.GroupOnGeom(v, k)

mesh.ExportUNV("meshCheckTestMesh.unv")		

MQC = MeshQualityCheck( mesh, 70, 2 )
MQC.checkMesh()
