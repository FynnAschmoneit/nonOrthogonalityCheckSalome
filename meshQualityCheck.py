# meshQualityCheck.py
from __future__ import division
import random
import re

import salome
salome.salome_init()
import GEOM, SMESH, SALOMEDS

from salome.geom import geomBuilder
geompy = geomBuilder.New(salome.myStudy)

from salome.smesh import smeshBuilder
smesh =  smeshBuilder.New(salome.myStudy)


class MeshQualityCheck:

	def __init__(self, MESH):
		self.mesh = MESH
		self.progressOutput = True
		self.nonOrthogonalityCheck = True
		self.skewnessCheck = True
		self.groupsForFailedVolumePairs = True
		
		self.outputPreSting = "MeshQualityCheck:\t\t"
		self.sharedFaceDict = {}
		self.nbInternalFaces = 0

		self.nonOrthogonalVolPairs = []
		self.nonOrthogonalThreshold = 70
		self.avNonOrth = 0
		self.maxNonOrth = 0

		self.tooSkewVolPairs = []
		self.skewnessThreshold = 2
		self.avSkew = 0
		self.maxSkew = 0

		random.seed(1)		# for coloring of volume pair groups

		if(self.progressOutput):
			print self.outputPreSting + "finding internal faces"

		volIDs = self.mesh.GetElementsByType( SMESH.VOLUME )
		for v in volIDs:
			nbF = self.mesh.ElemNbFaces( v )
			for f in range(0,nbF):
				vFNodes = self.mesh.GetElemFaceNodes( v, f )
				dictKey = "%s" % sorted(vFNodes)
				if dictKey not in self.sharedFaceDict:
					self.sharedFaceDict[ dictKey ] = [ v ]
				else:
					self.sharedFaceDict[ dictKey ].append( v )

		self.nbInternalFaces = len(self.sharedFaceDict) 

	def COMsLineDist(self, v_ar):
		bC1 = self.mesh.BaryCenter(v_ar[0])
		p1 = geompy.MakeVertex( bC1[0], bC1[1], bC1[2] )
		bC2 = self.mesh.BaryCenter(v_ar[1])
		p2 = geompy.MakeVertex( bC2[0], bC2[1], bC2[2] )
		line = geompy.MakeLineTwoPnt(p1, p2)
		centerDistance = geompy.MinDistance(p1, p2)
		return [line, centerDistance]

	def faceNormalDir(self, faceCoord):
		p1 = geompy.MakeVertex(faceCoord[0][0], faceCoord[0][1], faceCoord[0][2])
		p2 = geompy.MakeVertex(faceCoord[1][0], faceCoord[1][1], faceCoord[1][2])
		p3 = geompy.MakeVertex(faceCoord[2][0], faceCoord[2][1], faceCoord[2][2])
		v1 = geompy.MakeVector(p1, p2) 
		v2 = geompy.MakeVector(p2, p3)
		v3 = geompy.CrossProduct(v1, v2)
		return v3

	def checkNonOrth(self, COMsLine, fNorm, volPair):
		ang = geompy.GetAngle(COMsLine, fNorm)
		self.avNonOrth = self.avNonOrth + ang
		if ang > self.nonOrthogonalThreshold:
			self.nonOrthogonalVolPairs.append( volPair )
		if ang > self.maxNonOrth:
			self.maxNonOrth = ang

	def COMofFace(self, faceCoord, nbNodesOfFace):
		comFaceComp = [0]*3
		for i in range(0,3):
			for j in range(0, nbNodesOfFace):
				comFaceComp[i] = comFaceComp[i] + faceCoord[j][i]
		comFaceComp = [ comFaceComp[i]/nbNodesOfFace for i in range(0,3) ]
		return geompy.MakeVertex( comFaceComp[0], comFaceComp[1], comFaceComp[2] )
	
	def checkSkewness(self, pFCom, COMsLine, centerDistance, volPair):
		skewness = geompy.MinDistance(pFCom, COMsLine)/centerDistance
		self.avSkew = self.avSkew + skewness

		if skewness > self.skewnessThreshold:
			self.tooSkewVolPairs.append( volPair )
		if skewness > self.maxSkew:
			self.maxSkew = skewness

	def calcAverages(self):
		if( self.nonOrthogonalityCheck ):	
			self.avNonOrth = self.avNonOrth/len(self.sharedFaceDict)
		if( self.skewnessCheck ):
			self.avSkew = self.avSkew/len(self.sharedFaceDict)

	def printStats(self):
		print self.outputPreSting + "no faces: ", self.nbInternalFaces

		if( self.nonOrthogonalityCheck ):
			print self.outputPreSting + "non-orthogonality threshold: ", self.nonOrthogonalThreshold
			print self.outputPreSting + "no of non-orthogonal faces: ", len(self.nonOrthogonalVolPairs)
			print self.outputPreSting + "average non-orthogonality: ", self.avNonOrth
			print self.outputPreSting + "max non-orthogonality: ", self.maxNonOrth
		if( self.skewnessCheck ):
			print self.outputPreSting + "skewness threshold: ", self.skewnessThreshold
			print self.outputPreSting + "number of skew faces: ", len(self.tooSkewVolPairs)
			print self.outputPreSting + "average skewness: ", self.avSkew
			print self.outputPreSting + "max skewnes: ", self.maxSkew

	def createGroupsOfPairs(self):
		if( self.nonOrthogonalityCheck ):
			groupCounter = 0
			for v in self.nonOrthogonalVolPairs:
				interimGroup = self.mesh.GetMesh().CreateGroup(SMESH.VOLUME, "non-orthogonal pair "+ str(groupCounter) )
				color = SALOMEDS.Color( random.random(), random.random(), random.random() )
				interimGroup.Add( [ v[0], v[1] ] )
				interimGroup.SetColor(color)
				groupCounter = groupCounter + 1
		

		if( self.skewnessCheck ):
			groupCounter = 0
			for v in self.tooSkewVolPairs:
				interimGroup = self.mesh.GetMesh().CreateGroup(SMESH.VOLUME, "skew pair "+ str(groupCounter) )
				color = SALOMEDS.Color( random.random(), random.random(), random.random() )
				interimGroup.Add( [ v[0], v[1] ] )
				interimGroup.SetColor(color)
				groupCounter = groupCounter + 1
		
	def checkMesh(self):
		if(self.progressOutput):
			print self.outputPreSting + "checking mesh"

		faceIterNumber = 1
		for k, v in self.sharedFaceDict.iteritems():
			if len(v) == 2:	
			# k is the list of face node Ids for internal faces if corresponding volume vector consists of 2 elements 

				# calculation line connecting volumes' centers of mass and its distance
				COMsLine, centerDistance = self.COMsLineDist(v)
				
				# convert string of face node Ids k to list of corresp coordinates
				faceIDList = map(int, re.findall(r'\d+', k))
				nbNodesOfFace = len(faceIDList)
				faceCoord = [ self.mesh.GetNodeXYZ( faceIDList[i] ) for i in range(0, nbNodesOfFace) ]

				if( self.nonOrthogonalityCheck ):
				# calculating face normal direction
					fNorm = self.faceNormalDir(faceCoord)
				#calculating the angle between faceNormal and com-connecting line
					self.checkNonOrth(COMsLine, fNorm, v)

				if( self.skewnessCheck ):
				#calculation skewness as the minimum distance between the com-connecting line and com of the face
					pFCom = self.COMofFace(faceCoord, nbNodesOfFace)
					self.checkSkewness(pFCom, COMsLine, centerDistance, v)
			
			relProg = faceIterNumber/self.nbInternalFaces*100
			if( self.progressOutput and relProg % 5 < 99/self.nbInternalFaces  ):
				print self.outputPreSting + "progress: " +str(relProg//1) +"%"
			faceIterNumber = faceIterNumber + 1

		self.calcAverages()	
		self.printStats()
		if( self.groupsForFailedVolumePairs ):
			self.createGroupsOfPairs()
		
	

