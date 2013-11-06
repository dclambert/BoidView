################################################################################
#
#     boids.py
#
#     Python implementation of Craig Reynolds flocking alg.
#
#     Copyright (c) 2006, David C. Lambert
#
################################################################################

import sys
from math import *
from random import seed, uniform
from operator import add,sub

seed()
Range3 = range(3)

class Boid:
	_maxV = 10.0
	_boundingBox = ((-10.0, 10.0),) * 3

	def SetMaxVelocity(cls, v):
		Boid._maxV = v

	def MaxVelocity(cls):
		return Boid._maxV

	def SetBoundingBox(cls, b):
		Boid._boundingBox = b

	def BoundingBox(cls):
		return Boid._boundingBox

	MaxVelocity = classmethod(MaxVelocity)
	SetMaxVelocity = classmethod(SetMaxVelocity)
	BoundingBox = classmethod(BoundingBox)
	SetBoundingBox = classmethod(SetBoundingBox)

	def __str__(self):
		return str((self.pos, self.vel, self.accel))

	def __init__(self):
		self.oldTheta = None
		self.pos = [0.0, 0.0, 0.0]
		self.pos[0] = uniform(Boid._boundingBox[0][0], Boid._boundingBox[0][1])
		self.pos[1] = uniform(Boid._boundingBox[1][0], Boid._boundingBox[1][1])
		self.pos[2] = uniform(Boid._boundingBox[2][0], Boid._boundingBox[2][1])
		self.vel = [uniform(-5,5),uniform(-5,5),uniform(-5,5)]
		self.oldVel = self.vel
		self.magic = [uniform(0.8, 1.2),] * 3
		self.accel = [uniform(-5,5),uniform(-5,5),uniform(-5,5)]

	def clearAccelerations(self):
		self.accel = [0.0, 0.0, 0.0]

	def getDiffVectFrom(self, b):
		return [self.pos[i] - b.pos[i] for i in Range3]

	def applyMovements(self, mo):
		pX, pY, pZ = self.pos
		vX, vY, vZ = self.vel
		aX, aY, aZ = self.accel

		vX += aX * self.magic[0]
		vY += aY * self.magic[1]
		vZ += aZ * self.magic[2]

		vMag = 1.0e-06 + sqrt(vX*vX + vY*vY + vZ*vZ)
		if (vMag <= Boid._maxV):
			self.vel = [vX, vY, vZ]
		else:
			fact = Boid._maxV/vMag
			self.vel = [vX*fact, vY*fact, vZ*fact]

		self.vel = [(mo*self.oldVel[i] + (1.0-mo)*self.vel[i]) for i in Range3]

		pX = (pX + self.vel[0])
		pY = (pY + self.vel[1])
		pZ = (pZ + self.vel[2])

		minX, maxX = Boid._boundingBox[0]
		if (pX > maxX): pX = minX
		elif (pX < minX): pX = maxX

		minY, maxY = Boid._boundingBox[1]
		if (pY > maxY): pY = minY
		elif (pY < minY): pY = maxY

		minZ, maxZ = Boid._boundingBox[2]
		if (pZ > maxZ): pZ = minZ
		elif (pZ < minZ): pZ = maxZ

		self.pos = [pX, pY, pZ]
		self.oldVel = self.vel


class Flock:
	def __init__(self):
		self.maxV = 7.0
		self.minV = 1.0
		self.distExp = 2.2
		self.accLimit = 4
		self.avoidFact = 10.0
		self.matchFact = 0.7
		self.targetFact = 20.0
		self.centerFact = 0.1
		self.minRadius = 45.0
		self.minRadiusSq = pow(self.minRadius, self.distExp)
		self.momentum = 0.20
		self.minDist = 0.05
		self.distComp = 10.0   # adjust to boid diameter
		self.goalPos = [0.0, 0.0, 0.0]
		self.boidLambda = lambda x : Boid()

		self.nBoids = 0
		self.boids = []
		self.bbox = None
		self.mids = None

	def initBoids(self, n, bbox):
		self.bbox = bbox
		self.nBoids = n

		Boid.SetMaxVelocity(self.maxV)
		Boid.SetBoundingBox(self.bbox)

		self.ranges = map(lambda x: x[1]-x[0], self.bbox)
		self.mids = map(lambda x: x[0]+(x[1]-x[0])/2.0, self.bbox)
		self.boids = map(self.boidLambda, xrange(n))

		self.newGoal()

	def newGoal(self):
		for i in Range3:
			self.goalPos[i] = uniform(self.mids[i]-self.ranges[i]*.4, self.mids[i]+self.ranges[i]*.4)

	def applyMovements(self):
		for b in self.boids:
			b.applyMovements(self.momentum)

	def calcDists(self, posDiff):
		dist = sqrt(reduce(add, [x*x for x in posDiff])) - self.distComp
		if (dist < 0.0): dist = 0.1
		adjDist = pow(dist, self.distExp)
		return (dist, adjDist)

	def correctDiffForBBox(self, posDiff):
		for i in Range3:
			if (posDiff[i] > 0.9 * self.ranges[i]):
				posDiff[i] = self.ranges[i] - posDiff[i]

	def computeGroupVectors(self, refBoid):
		refBoidId = id(refBoid)
		neighborCount = 0
		avgVel = [0.0, 0.0, 0.0]
		centroid = [0.0, 0.0, 0.0]
		avoidance = [0.0, 0.0, 0.0]
		for j in range(self.nBoids):
			testBoid = self.boids[j]

			if (id(testBoid) == refBoidId): continue

			posDiff = [refBoid.pos[i] - testBoid.pos[i] for i in Range3]
			#self.correctDiffForBBox(posDiff);

			dist, adjDist = self.calcDists(posDiff)
			if (adjDist > self.minRadiusSq): continue

			neighborCount += 1

			avgVel = [avgVel[i] + testBoid.vel[i] for i in Range3]
			centroid = [centroid[i] + testBoid.pos[i] for i in Range3]

			avoidance = [avoidance[i] + posDiff[i]/adjDist for i in Range3]


		if 0:
			for i in Range3:
				posDiff[i] = refBoid.pos[i] - self.bbox[i][0]
			dist, adjDist = self.calcDists(posDiff)
			if (1 or adjDist <= self.minRadiusSq):
				avoidance = [avoidance[i] + 50*posDiff[i]/adjDist for i in Range3]

			for i in Range3:
				posDiff[i] = refBoid.pos[i] - self.bbox[i][1]
			dist, adjDist = self.calcDists(posDiff)
			if (1 or adjDist <= self.minRadiusSq):
				avoidance = [avoidance[i] + 50*posDiff[i]/adjDist for i in Range3]

		if (neighborCount > 0):
			avgVel = [x/neighborCount for x in avgVel]
			centroid = [x/neighborCount for x in centroid]

		return (neighborCount, avoidance, avgVel, centroid)

	def computeAccelerations(self):
		# avoid flockmates
		for i in range(self.nBoids):
			neighborCount = 0
			avgVel = [0.0, 0.0, 0.0]
			centroid = [0.0, 0.0, 0.0]
			avoidance = [0.0, 0.0, 0.0]
			refBoid = self.boids[i]
			refBoid.clearAccelerations()

			(neighborCount, avoidance, avgVel, centroid) = self.computeGroupVectors(refBoid)

			# flockmate avoidance - CHECK
			avoidAccel = [x*self.avoidFact for x in avoidance]
			refBoid.accel = [refBoid.accel[i] + avoidAccel[i] for i in Range3]

			accMag = sqrt(reduce(add, [x*x for x in refBoid.accel]))
			if (neighborCount > 0 and accMag < self.accLimit):
				# velocity matching - CHECK
				matchAccel = [(avgVel[i] - refBoid.vel[i]) * self.matchFact for i in Range3]
				refBoid.accel = [refBoid.accel[i] + matchAccel[i] for i in Range3]

				# flock centering - CHECK
				accMag = sqrt(reduce(add, [x*x for x in refBoid.accel]))
				if (accMag < self.accLimit):
					centerAccel = [(centroid[i] - refBoid.pos[i])*self.centerFact for i in Range3]
					refBoid.accel = [refBoid.accel[i] + centerAccel[i] for i in Range3]

			# target following - CHECK
			accMag = sqrt(reduce(add, [x*x for x in refBoid.accel]))
			if (accMag < self.accLimit):
				posDiff = [self.goalPos[i] - refBoid.pos[i] for i in Range3]
				#self.correctDiffForBBox(posDiff);

				dist, adjDist  = self.calcDists(posDiff)

				if (dist > self.minRadius):
					targetAccel = [x*self.targetFact/adjDist for x in posDiff]
					refBoid.accel = [refBoid.accel[i] + targetAccel[i] for i in Range3]

