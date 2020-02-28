



class Media:

	def __init__(self, volume, carrying_capacity=2.0e9):

		self.volume				= volume
		self.carryingCapacity 	= carrying_capacity


##################################################
##################################################


class Plasmid:

	def __init__(self, name, copy_number=1, marker=''):

		self.name 				= name

		self.copyNumber 		= copy_number
		self.marker 			= marker

##################################################
##################################################


class Strain:

	def __init__(self, name, max_growth_rate, optimal_temp=37.0, mean_lag_exit_time=1.5, stdev_lag_exit_time=0.125, marker='', plasmids=[]):

		self.name 				= name

		self.maxGrowthRate 		= max_growth_rate
		self.optimalTemp		= optimal_temp
		self.meanLagExitTime 	= mean_lag_exit_time
		self.stdevLagExitTime 	= stdev_lag_exit_time

		self.marker 			= marker

		self.plasmids 			= plasmids


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