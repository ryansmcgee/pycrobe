from __future__ import division
from decimal import Decimal

import random
import numpy as numpy
import pandas as pandas
import scipy.stats


##################################################
##################################################


class Media:

	def __init__(self, volume, carrying_capacity=2.0e9):

		self.volume				= volume
		self.carryingCapacity 	= carrying_capacity


##################################################
##################################################


class Strain:

	def __init__(self, name, max_growth_rate, optimal_temp=37.0, mean_lag_exit_time=1.5, stdev_lag_exit_time=0.125, marker=''):

		self.name 				= name

		self.maxGrowthRate 		= max_growth_rate
		self.optimalTemp		= optimal_temp
		self.meanLagExitTime 	= mean_lag_exit_time
		self.stdevLagExitTime 	= stdev_lag_exit_time

		self.marker 			= marker


##################################################
##################################################


class Inoculum:

	def __init__(self, strain, cell_count, growth_phase='stationary', growth_cycle_timer=0):

		self.strain 			= strain
		self.cellCount 			= cell_count

		self.growthPhase 		= growth_phase # 'stationary', 'lag', 'exponential'
		self.growthCycleTimer	= growth_cycle_timer


	def growthRate(self, temp, saturationCoeff, dt):
		r 	= self.strain.maxGrowthRate 	# intrinsic max exponential growth rate
		r 	*= temp/self.strain.optimalTemp	# growth rate dependence on temperature
		r 	*= scipy.stats.norm.cdf(self.growthCycleTimer, loc=self.strain.meanLagExitTime, scale=self.strain.stdevLagExitTime) 	# modeling lag phase exit (cells exit lag phase (r=0) and enter exponential phase (r=maxrate*otherfactors) at lag-exit-times that are normally distributed from cell to cell; this is approximated by the population's average growth rate following the normal distn cdf centered at the mean lag exit time
		r 	*= saturationCoeff 	# growth rate dependence on carrying capacity
		r 	= numpy.random.normal(r, r*0.02)	# additional random noise

		self.growthCycleTimer += dt
		if(saturationCoeff <= 0.01):
			self.growthPhase = 'stationary'
		elif(self.growthPhase == 'stationary' and saturationCoeff > 0.01):
			self.growthPhase = 'lag'
			self.growthCycleTimer = 0.0
		elif(self.growthPhase == 'lag' and (self.strain.maxGrowthRate-r)/self.strain.maxGrowthRate < 0.01):
			self.growthPhase = 'exponential'

		# print str(self.strain.name) + ' ' + str(self.strain.maxGrowthRate) + ' ' + str(r)+ ' ' + str(self.growthPhase) + ' ' + str(self.growthCycleTimer)
		
		return r


	def incubate(self, temp, dt, saturationCoeff):

		# No drug case
		dN 	= self.growthRate(temp, saturationCoeff, dt) * self.cellCount
		self.cellCount += dN*dt


##################################################
##################################################


class Drug:

	def __init__(self, name, mg, spontaneous_degradation_rate=0):

		self.name 							= name
		self.mg 							= mg		
		self.spontaneous_degradation_rate 	= spontaneous_degradation_rate

	def incubate(self, temp, dt):
		# Do drug degradation
		return


##################################################
##################################################


class CultureVolume:

	def __init__(self, media=[], inoculums=[], drugs=[], name=""):

		self.name 		= name

		self.media	= []	# class instance list has to be declared this way to avoid list being shared among instances when mutable changes are made to one instance list, apparently
		for media in (media if isinstance(media, list) else [media]):
			self.media.append(media)

		self.inoculums	= []	# class instance list has to be declared this way to avoid list being shared among instances when mutable changes are made to one instance list, apparently
		for inoculum in (inoculums if isinstance(inoculums, list) else [inoculums]):
			self.inoculums.append(inoculum)

		self.drugs	= []	# class instance list has to be declared this way to avoid list being shared among instances when mutable changes are made to one instance list, apparently
		for drug in (drugs if isinstance(drugs, list) else [drugs]):
			self.drugs.append(drug)


	def totalVolume(self):
		return numpy.sum([m.volume for m in self.media]) if len(self.media)>0 else 0.0


	def totalCarryingCapacity(self):
		return numpy.sum([m.carryingCapacity*(m.volume/self.totalVolume()) for m in self.media]) if len(self.media)>0 else 0.0


	def getCellCounts(self):
		cellCounts = {}
		for inoculum in self.inoculums:
			cellCounts[inoculum.strain.name] = inoculum.cellCount
		return cellCounts


	def getCellDensities(self):
		cellDensities = {}
		for inoculum in self.inoculums:
			cellDensities[inoculum.strain.name] = inoculum.cellCount / self.totalVolume()
		return cellDensities


	def getInoculumFrequencies(self):
		frequencies = {}
		for inoculum, count in self.getCellCounts().iteritems():
			frequencies[inoculum] = count / self.totalCellCount()
		return frequencies


	def getMarkerCounts(self):
		markerCounts = {}
		for inoculum in self.inoculums:
			markerCounts[inoculum.strain.marker] = markerCounts.get(inoculum.strain.marker, 0) + inoculum.cellCount
		return markerCounts


	def getMarkerDensities(self):
		markerDensities = {}
		markerCounts = self.getMarkerCounts()
		# print markerCounts
		# exit()
		for marker, markerCount in markerCounts.iteritems():
			markerDensities[marker] = markerCount / self.totalVolume()
		return markerDensities


	def getDrugConcentrations(self):
		drugConcentrations = {}
		for drug in self.drugs:
			drugConcentrations[drug.name] = drug.mg / self.totalVolume()
		return drugConcentrations


	def totalCellCount(self):
		return numpy.sum([i.cellCount for i in self.inoculums])


	def totalCellDensity(self):
		return numpy.sum([i.cellCount/self.totalVolume() for i in self.inoculums])


	def add(self, addedCultureVolume):
		for addedMedia in addedCultureVolume.media:
			self.media.append(addedMedia)
		for addedInoculum in addedCultureVolume.inoculums:
			inoculumExistsInCulture = False
			for i, inoculum in enumerate(self.inoculums):
				if(inoculum.strain.name == addedInoculum.strain.name):
					self.inoculums[i].cellCount += addedInoculum.cellCount
					inoculumExistsInCulture = True
					break
			if(not inoculumExistsInCulture):
				self.inoculums.append(addedInoculum)
		for addedDrug in addedCultureVolume.drugs:
			drugExistsInCulture = False
			for i, drug in enumerate(self.drugs):
				if(drug.name == addedDrug.name):
					self.drugs[i].mg += addedDrug.mg
					drugExistsInCulture = True
					break
			if(not drugExistsInCulture):
				self.drugs.append(addedDrug)	


	def sample(self, volume, randomError=True, errorTolerancePct=2.0):
		presampleVolume = self.totalVolume()
		sampledVolume 	= numpy.random.normal(volume, volume*(errorTolerancePct/100.0))
		sampledMedia 	= []
		for media in self.media:
			sampledMediaVolume 	= media.volume * sampledVolume/presampleVolume
			media.volume 		-= sampledMediaVolume
			sampledMedia.append( Media(volume=sampledMediaVolume, carrying_capacity=media.carryingCapacity) )
		sampledInoculums = []
		for inoculum in self.inoculums:
			try:
				sampledCellCount 	= numpy.random.poisson( inoculum.cellCount * sampledVolume/presampleVolume )	# inoculum.cellCount * sampledVolume/presampleVolume )	# num cells sampled = draw from poisson distribution given expected num cells per sample volume
			except ValueError:
				sampledCellCount 	= 0
			inoculum.cellCount 		-= sampledCellCount
			sampledInoculums.append( Inoculum(strain=inoculum.strain, cell_count=sampledCellCount, growth_phase=inoculum.growthPhase, growth_cycle_timer=inoculum.growthCycleTimer) )
		sampledDrugs 	= []
		for drug in self.drugs:
			try:
				sampledDrugNg 	= numpy.random.poisson( drug.mg*1000000 * sampledVolume/presampleVolume )	# drug.cellCount*1000000 * sampledVolume/presampleVolume )	# num *ng* sampled = draw from poisson distribution given expected num ng's per sample volume
				sampledDrugMg	= sampledDrugNg/1000000
			except ValueError:
				sampledDrugMg 	= 0
			drug.mg 			-= sampledDrugMg
			sampledDrugs.append( Drug(name=drug.name, mg=sampledDrugMg) )
		# self.totalVolume() 		-= sampledVolume
		sampledCulture 		= CultureVolume(media=sampledMedia, inoculums=sampledInoculums, drugs=sampledDrugs, name=self.name+" sample")
		return sampledCulture


	def incubate(self, temp, dt):
		saturationCoeff = (self.totalCarryingCapacity() - self.totalCellDensity())/self.totalCarryingCapacity()
		for inoculum in self.inoculums:
			inoculum.incubate(temp=temp, dt=dt, saturationCoeff=saturationCoeff)


	def info(self):
		print "CultureVolume " + self.name + ":"
		print "--------------------------------"
		print "\tVolume\t\t\t= " + ("%.3f" % self.totalVolume()) + " mL"
		print "\tCarrying Capacity\t= " + ( '%.2E' % Decimal(str( self.totalCarryingCapacity() )) ) + " cfu/mL"
		print "\tDensity\t\t\t= " + ( '%.2E' % Decimal(str( self.totalCellDensity() )) ) + " cfu/mL"
		print "\tTotal Cell Count\t= " + ( '%.2E' % Decimal(str( self.totalCellCount() )) ) + " cfu"
		if(len(self.inoculums)>0):
			for inoculum in self.inoculums:
				print "\tInoculum " + inoculum.strain.name + ":\n" + "\t\tdensity\t\t= " + ( '%.2E' % Decimal(str( inoculum.cellCount/self.totalVolume() )) )+" cfu/mL" + "\n\t\tcell count\t= " + ( '%.2E' % Decimal(str( inoculum.cellCount )) )+" cfu"
		else:
			print "\t(no inoculums)"
		if(len(self.drugs)>0):
			for drug in self.drugs:
				print "\tDrug " + drug.name + ":\n" + "\t\tconcentration\t= " + ("%.5f" % self.getDrugConcentrations()[drug.name]) + " mg/mL"
		else:
				print "\t(no drugs)"


##################################################
##################################################


class DilutionSeries:

	def __init__(self, initial_culture, fold_dilution=10, num_serial_dilutions=10, transfer_volume=0.030):

		self.initialCulture 	= initial_culture
		self.foldDilution 		= fold_dilution
		self.numSerialDilutions = num_serial_dilutions
		self.transferVolume 	= transfer_volume
		
		self.dilutions 	 		= [None for d in range(self.numSerialDilutions+1)]
		self.dilutions[0] 		= self.initialCulture

		for d in range(1, self.numSerialDilutions+1):
			self.dilutions[d] = CultureVolume( 	media=[Media(volume=(self.transferVolume*self.foldDilution)-self.transferVolume, carrying_capacity=0)],
												name=self.initialCulture.name+"Serial Dilution -" + str(d)
											 )
			self.dilutions[d].add(self.dilutions[d-1].sample(self.transferVolume))

			
	def getDilution(self, dilution):
		d = numpy.abs(dilution)
		return self.dilutions[d]


	def getDensestCountableDilution(self):
		for d in range(self.numSerialDilutions):
			if(self.dilutions[d].totalCellCount()<400):
				return {'culture': self.dilutions[d], 'df': d}



##################################################
##################################################


class AgarPlate:

	def __init__(self, drugs=[], num_colony_landing_spots=10000, name=""):

		self.name 		= name
		self.inoculums	= []
		self.drugs 		= drugs 	# TODO: handle drug concs and effect on colonies forming baesd on cell MICs

		self.numColonyLandingSpots 	= num_colony_landing_spots
		self.colonyLandingSpots 	= [[] for n in range(num_colony_landing_spots)]


	def inoculate(self, addedCultureVolume):
		for addedInoculum in addedCultureVolume.inoculums:
			existingInoculum = next((i for i in self.inoculums if i.strain == addedInoculum.strain.name), None)
			if(existingInoculum):
				existingInoculum.cellCount += addedInoculum.cellCount
			else:
				self.inoculums.append(addedInoculum)

		for inoculum in self.inoculums:
			for cell in range(inoculum.cellCount):
				landingSpotIndex = numpy.random.randint(low=0, high=self.numColonyLandingSpots)
				if(not self.colonyLandingSpots[landingSpotIndex]):
					self.colonyLandingSpots[ landingSpotIndex ] = [inoculum.strain.marker]
				else:
					self.colonyLandingSpots[ landingSpotIndex ].append(inoculum.strain.marker)

		for l, landingSpotMarkers in enumerate(self.colonyLandingSpots):
				if(len(landingSpotMarkers) >= 1):
					self.colonyLandingSpots[l] = random.choice(landingSpotMarkers) # Of all the cells that landed on the same spot, one random one will become a visible colony
				else:
					self.colonyLandingSpots[l] = None


	def incubate(self, temp=37.0, dt=1):
		if(temp > 30):
			for l, landingSpotMarkers in enumerate(self.colonyLandingSpots):
				if(landingSpotMarkers != 0 and landingSpotMarkers != ''):
					self.colonyLandingSpots[l] = random.choice(landingSpotMarkers) # Of all the cells that landed on the same spot, one random one will become a visible colony


	def getCellCounts(self):
		cellCounts = {}
		for inoculum in self.inoculums:
			cellCounts[inoculum.strain.name] = inoculum.cellCount


	def totalCellCount(self):
		return numpy.sum([i.cellCount for i in self.inoculums])


	def getColonyCounts(self):
		colonyCounts = {}
		for marker in numpy.unique([inoculum.strain.marker for inoculum in self.inoculums]):
			colonyCounts[marker] = sum((L!=None and marker in L) for L in self.colonyLandingSpots)
		return colonyCounts


	def totalColonyCount(self):
		return sum(L!=None for L in self.colonyLandingSpots)


	def info(self):
		print "AgarPlate " + self.name + ":"
		print "--------------------------------"
		print "\tTotal Colonies = " + str(self.totalColonyCount())
		for marker in numpy.unique([inoculum.strain.marker for inoculum in self.inoculums]):
			print "\tMarker " + marker + ":\n" + "\t\tcell count\t= " + str( self.getColonyCounts()[marker] )
		# for inoculum in self.inoculums:
		# 	print "\tInoculum " + inoculum.strain.name + ":\n" + "\t\tcell count\t= " + str( self.getColonyCounts()[inoculum.strain.marker] )


##################################################
##################################################


class Incubator:

	def __init__(self, set_temp, temp_std_batch=0, temp_std_location=0, temp_std_transient=0):

		self.setTemp 			= set_temp
		self.tempSTD_batch		= temp_std_batch
		self.tempSTD_location 	= temp_std_location
		self.tempSTD_transient 	= temp_std_transient

	def incubate(self, cultures, time, dt, record_growth_curves=False):

		if(not isinstance(cultures,(list, numpy.ndarray))):
			cultures = [cultures]

		baseTemp_batch 		= numpy.random.normal(self.setTemp, self.tempSTD_batch) 

		baseTemp_locations 	= [ numpy.random.normal(baseTemp_batch, self.tempSTD_location) for culture in cultures ]

		effectiveTemps 		= baseTemp_locations

		growthCurves 	= []

		t=0
		for t in numpy.arange(start=0, stop=time, step=dt):
			
			effectiveTemps 	+= numpy.random.normal(0, self.tempSTD_transient, len(effectiveTemps))
			
			for c, culture in enumerate(cultures):
				culture.incubate(temp=effectiveTemps[c], dt=dt)

				if(record_growth_curves):
					growthCurves.append( {'culture':culture.name, 'time':t+dt, 'temp':effectiveTemps[c], 'density':culture.totalCellDensity()} )

		if(record_growth_curves):
			return pandas.DataFrame(growthCurves)




##################################################
##################################################


class FlowCytometer:

	def __init__(self, name=""):

		self.name 		= name
		
		# self.eventSlotsPerUl 	= event_slots_per_ul
		# self.eventSlots 	= [[] for n in range(num_colony_landing_spots)]


	def read(self, sample, read_volume, event_slots_per_ul=10000):
		numEventSlots 	= event_slots_per_ul * int(read_volume*1000)
		eventSlots 		= [[] for n in range(numEventSlots)]
		for inoculum in sample.inoculums:
			for cell in range(inoculum.cellCount):
				eventSlotIndex = numpy.random.randint(low=0, high=numEventSlots)
				if(not eventSlots[eventSlotIndex]):
					eventSlots[ eventSlotIndex ] = [inoculum.strain.marker]
				else:
					eventSlots[ eventSlotIndex ].append(inoculum.strain.marker)
		cellEventCounts = {}
		for eventMarkers in eventSlots:
			overallEventMarker = '-'.join(numpy.unique(eventMarkers))
			if overallEventMarker in cellEventCounts:
				cellEventCounts[overallEventMarker] += 1
			else:
				cellEventCounts[overallEventMarker] = 1
		try:
			del cellEventCounts['']
		except KeyError:
			pass
		return cellEventCounts


		

		




















