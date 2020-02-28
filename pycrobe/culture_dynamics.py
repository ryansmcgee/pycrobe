from __future__ import division

import numpy
import pandas
import scipy.stats

# from standard import *
# from betalactamase import *


##################################################
##################################################


class CultureDynamics(object):

	def __init__(self):
		
		self.data 		= pandas.DataFrame(columns=['run', 'strain', 't', 't_run', 'N', 'r', 'S'])

		self.num_runs 	= 0


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	def run(self, time, dt, inoculums, media, temp, noise_r=0.02, downsample_output_dt=None):
		
		I 	= len(inoculums)

		T 	= int(time/dt)

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# LOAD PARAMETER VALUES: 
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		rho 			= numpy.zeros(I)
		k_r 			= numpy.zeros(I)
		sigma 			= numpy.zeros(I)

		optimalTemp			= numpy.zeros(I)
		meanLagExitTime 	= numpy.zeros(I)
		stdevLagExitTime 	= numpy.zeros(I)

		#--------------------------------------------------

		for i, inoculum in enumerate(inoculums):

			rho[i] 			= inoculum.strain.maxGrowthRate
			k_r[i] 			= inoculum.strain.halfmaxGrowthNutrientConc
			sigma[i] 		= inoculum.strain.nutrientConsumptionRate

			optimalTemp[i] 		= inoculum.strain.optimalTemp
			meanLagExitTime[i] 	= inoculum.strain.meanLagExitTime
			stdevLagExitTime[i] = inoculum.strain.stdevLagExitTime


		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# LOAD INITIAL VARIABLE VALUES: 
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		N 		= numpy.zeros(shape=(T+1,I))	
		r 		= numpy.zeros(shape=(T+1,I))	

		growthPhase 		= [''] * I
		growthCycleTimer	= numpy.zeros(I)

		S 		= numpy.zeros(shape=(T+1,1))		

		#--------------------------------------------------

		for i, inoculum in enumerate(inoculums):

			N[0][i] 	= inoculum.cellCount
			r[0][i] 	= inoculum.growthRate

			growthPhase[i] 		= inoculum.growthPhase
			growthCycleTimer[i] = inoculum.growthCycleTimer

		#--------------------------------------------------

		S[0] = media.nutrient.concentration * media.volume


		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# RUN GROWTH MODEL:
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		for t in range(T):

			dN 			= r[t]*N[t]
			N[t+1] 		= N[t] + dN*dt

			dS 			= -sigma * numpy.sum(r[t]*N[t])
			S[t+1] 		= max(S[t] + dS*dt, 0)

			# Modulate the effective max growth rate by growth phase and temperature:
			r_max 		= numpy.zeros(I)
			for i in range(I):
				r_max[i] = rho[i] 					# intrinsic max exponential growth rate
				r_max[i] *= scipy.stats.norm.cdf(growthCycleTimer[i], loc=meanLagExitTime[i], scale=stdevLagExitTime[i]) 	# modeling lag phase exit (cells exit lag phase (r=0) and enter exponential phase (r=maxrate*otherfactors) at lag-exit-times that are normally distributed from cell to cell; this is approximated by the population's average growth rate following the normal distn cdf centered at the mean lag exit time
				r_max[i] *= temp/optimalTemp[i]		# growth rate dependence on temperature
				# Apply additional random noise:
				r_max[i] = numpy.random.normal(r_max[i], r_max[i]*noise_r)
			# Calculate the current growth rate as a function of nutrient availability:
			r[t+1] 		= r_max * ( S[t] / (k_r + S[t]) )
			# Update the current growth cycle phase:
			for i in range(I):
				growthCycleTimer[i] += dt
				if( (S[t] / (k_r + S[t])) <= 0.01 ):
					growthPhase[i] = 'stationary'
				elif( growthPhase[i] == 'stationary' and (S[t] / (k_r + S[t])) > 0.01 ):
					growthPhase[i] = 'lag'
					growthCycleTimer[i] = 0.0
				elif( growthPhase[i] == 'lag' and (rho[i]-r[i])/rho[i] < 0.5 ):
					growthPhase[i] = 'exponential'


		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# SAVE OUT FINAL VARIABLE VALUES: 
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		for i, inoculum in enumerate(inoculums):

			inoculum.cellCount 	= N[-1][i]
			inoculum.growthRate = r[-1][i]

			inoculum.growthPhase 		= growthPhase[i]
			inoculum.growthCycleTimer	= growthCycleTimer[i]
			
		media.nutrient.concentration = S[-1] / media.volume


		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# UPDATE THE DATAFRAME HOLDING THE VARIABLE TIME SERIES:
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		for i, strainName in enumerate([inoculum.strain.name for inoculum in inoculums]):
			t_vals 			= numpy.arange(start=0, stop=time+dt, step=dt)
			t_vals_macro 	= t_vals if self.num_runs==0 else t_vals+self.data['t'].max() 	# in terms of an overall timer, assume that each run starts immediately where the previous run left off
			run_timeseries 	= {
								'run': [self.num_runs]*(T+1),
								'strain': [strainName]*(T+1),
								't': t_vals_macro,
								't_run': t_vals,
								'N': N[:,i],
								'r': r[:,i],
								'S': S[:,i]
							  }

			strain_run_data = pandas.DataFrame.from_dict(run_timeseries)

			self.data = self.data.append(strain_run_data)
			self.data = self.data.drop_duplicates(subset=['t', 'strain'], keep='first')


		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		self.num_runs 	+= 1

		return


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	def figure(self):
		import matplotlib.pyplot as pyplot
		import seaborn

		fig = pyplot.figure(figsize=(16,9))
		ax1 = pyplot.subplot(221)
		ax2 = pyplot.subplot(223)
		ax3 = pyplot.subplot(222)

		seaborn.lineplot(data=self.data, ax=ax1, x="t", y="N", hue="strain")

		seaborn.lineplot(data=self.data, ax=ax2, x="t", y="r", hue="strain")

		seaborn.lineplot(data=self.data, ax=ax3, x="t", y="S")

		seaborn.set_style("ticks")
		seaborn.despine()

		fig.tight_layout()

		pyplot.show()


###################################################
###################################################


class BetaLactamaseDynamics(CultureDynamics):

	def __init__(self):
		super(BetaLactamaseDynamics, self).__init__()
		return


	def run(self, time, dt, inoculums, media, temp, noise_r=0.02, downsample_output_dt=None):

		I 	= len(inoculums)

		T 	= int(time/dt)

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# LOAD PARAMETER VALUES: 
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		rho 			= numpy.zeros(I)
		k_r 			= numpy.zeros(I)
		lamda 			= numpy.zeros(I)
		k_l 			= numpy.zeros(I)
		eta 			= numpy.zeros(I)
		alpha_B 		= numpy.zeros(I)
		delta_B 		= numpy.zeros(I)
		delta_Bext 		= numpy.zeros(I)
		epsilon_B 		= numpy.zeros(I)
		xi_B 			= numpy.zeros(I)
		delta_A 		= numpy.zeros(I)
		epsilon_A 		= numpy.zeros(I)
		xi_A 			= numpy.zeros(I)
		Vmax_B 			= numpy.zeros(I)
		k_B 			= numpy.zeros(I)
		Vmax_Bext 		= numpy.zeros(I)
		k_Bext 			= numpy.zeros(I)
		Vol 			= numpy.zeros(I)
		sigma 			= numpy.zeros(I)

		optimalTemp			= numpy.zeros(I)
		meanLagExitTime 	= numpy.zeros(I)
		stdevLagExitTime 	= numpy.zeros(I)

		#--------------------------------------------------

		for i, inoculum in enumerate(inoculums):

			rho[i] 			= inoculum.strain.maxGrowthRate
			k_r[i] 			= inoculum.strain.halfmaxGrowthNutrientConc
			lamda[i] 		= inoculum.strain.maxLysisRate
			k_l[i] 	 		= inoculum.strain.halfmaxLysisDrugConc
			eta[i] 			= inoculum.strain.lysisHillCoefficient
			alpha_B[i] 		= inoculum.strain.betalactamaseProductionRate
			delta_B[i] 		= inoculum.strain.betalactamase.decayRate_intracellular
			delta_Bext[i]	= inoculum.strain.betalactamase.decayRate_extracellular
			epsilon_B[i] 	= inoculum.strain.betalactamaseLeakRate
			xi_B[i] 		= inoculum.strain.betalactamaseDebrisSinkFraction
			epsilon_A[i] 	= inoculum.strain.betalactamDiffusionRate
			xi_A[i] 		= inoculum.strain.betalactamDebrisSinkFraction
			Vmax_B[i] 		= inoculum.strain.betalactamase.maxHydrolysisRate
			k_B[i] 			= inoculum.strain.betalactamase.halfmaxHydrolysisDrugConc
			Vmax_Bext[i] 	= inoculum.strain.betalactamase.maxHydrolysisRate
			k_Bext[i] 		= inoculum.strain.betalactamase.halfmaxHydrolysisDrugConc
			Vol[i] 			= inoculum.strain.periplasmVolume
			sigma[i] 		= inoculum.strain.nutrientConsumptionRate

			optimalTemp[i] 		= inoculum.strain.optimalTemp
			meanLagExitTime[i] 	= inoculum.strain.meanLagExitTime
			stdevLagExitTime[i] = inoculum.strain.stdevLagExitTime

		#--------------------------------------------------

		for drug in media.drugs:
			if(isinstance(drug, BetaLactam)):
				delta_A = drug.decayRate
				# Assume Media objects hold only one drug object per type (BetaLactam)
				break


		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# LOAD INITIAL VARIABLE VALUES: 
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		N 		= numpy.zeros(shape=(T+1,I))	
		r 		= numpy.zeros(shape=(T+1,I))	
		l 		= numpy.zeros(shape=(T+1,I))	
		B 		= numpy.zeros(shape=(T+1,I))	
		Bext 	= numpy.zeros(shape=(T+1,I))	
		d_B 	= numpy.zeros(shape=(T+1,I))	
		d_Bext 	= numpy.zeros(shape=(T+1,I))	
		A 		= numpy.zeros(shape=(T+1,I))	
		Aext 	= numpy.zeros(shape=(T+1,1)) 	
		S 		= numpy.zeros(shape=(T+1,1))		

		growthPhase 		= [''] * I
		growthCycleTimer	= numpy.zeros(I)

		#--------------------------------------------------

		for i, inoculum in enumerate(inoculums):

			N[0][i] 	= inoculum.cellCount
			r[0][i] 	= inoculum.growthRate
			l[0][i] 	= inoculum.lysisRate

			growthPhase[i] 		= inoculum.growthPhase
			growthCycleTimer[i] = inoculum.growthCycleTimer

			for solute in inoculum.periplasm.solutes:
				if(isinstance(solute, BetaLactamase)):
					B[0][i] 	= solute.concentration
					d_B[0][i] 	= solute.hydrolysisRate

					# Find the corresponding extracellular BetaLactamase:
					matchingBlaFound = False
					for solute_ext in media.solutes:
						if(isinstance(solute_ext, BetaLactamase) and solute == solute_ext):
							Bext[0][i] 	= solute_ext.concentration
							d_Bext[0][i] 	= solute_ext.hydrolysisRate
							matchingBlaFound = True
							break
					if(not matchingBlaFound): 
						# Make new object for Bla type that doesn't currently exist in extracellular media:
						betalactamase_ext = BetaLactamase(name=solute.name)
						betalactamase_ext.__dict__.update(solute.__dict__)
						betalactamase_ext.concentration = 0.0
						media.solutes.append(betalactamase_ext)

					# Assume Media objects hold only one solute object per type (BetaLactamase)
					break

			A[0][i] = 0.0
			for drug in inoculum.periplasm.drugs:
				if(isinstance(drug, BetaLactam)):
					A[0][i] 	= drug.concentration
					# Assume Media objects hold only one drug object per type (BetaLactam)
					break

		#--------------------------------------------------

		Aext[0] = 0.0
		for drug in media.drugs:
			if(isinstance(drug, BetaLactam)):
				Aext[0] = drug.concentration
				# Assume Media objects hold only one drug object per type (BetaLactam)
				break

		S[0] = media.nutrient.concentration * media.volume


		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# RUN GROWTH MODEL:
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		for t in range(T):

			dN 			= (r[t] - l[t])*N[t]
			N[t+1] 		= N[t] + dN*dt

			dB 			= -delta_B*B[t] - r[t]*B[t] - epsilon_B*(B[t] - Bext[t]) + alpha_B
			B[t+1] 		= B[t] + dB*dt

			dBext 		= -delta_Bext*Bext[t] + numpy.sum(epsilon_B*(B[t] - Bext[t])) + l[t]*xi_B*Vol*B[t]*N[t]
			Bext[t+1] 	= Bext[t] + dBext*dt

			dA 			= -delta_A*A[t] + epsilon_A*(Aext[t] - A[t]) - d_B[t]*B[t]
			A[t+1] 		= A[t] + dA*dt

			dAext 		= -delta_A*Aext[t] - numpy.sum(epsilon_A*(Aext[t] - A[t])) - numpy.sum(d_Bext*Bext[t]) + l[t]*xi_A*Vol*A[t]*N[t]
			Aext[t+1] 	= Aext[t] + dAext*dt

			dS 			= -sigma * numpy.sum(r[t]*N[t])
			S[t+1] 		= S[t] + dS*dt

			# Modulate the effective max growth rate by growth phase and temperature:
			r_max 		= numpy.zeros(I)
			for i in range(I):
				r_max[i] = rho[i] 					# intrinsic max exponential growth rate
				r_max[i] *= scipy.stats.norm.cdf(growthCycleTimer[i], loc=meanLagExitTime[i], scale=stdevLagExitTime[i]) 	# modeling lag phase exit (cells exit lag phase (r=0) and enter exponential phase (r=maxrate*otherfactors) at lag-exit-times that are normally distributed from cell to cell; this is approximated by the population's average growth rate following the normal distn cdf centered at the mean lag exit time
				r_max[i] *= temp/optimalTemp[i]		# growth rate dependence on temperature
				# Apply additional random noise:
				r_max[i] = numpy.random.normal(r_max[i], r_max[i]*noise_r)
			# Calculate the current growth rate as a function of nutrient availability:
			r[t+1] 		= r_max * ( S[t] / (k_r + S[t]) )
			# Update the current growth cycle phase:
			for i in range(I):
				growthCycleTimer[i] += dt
				if( (S[t] / (k_r + S[t])) <= 0.01 ):
					growthPhase[i] = 'stationary'
				elif( growthPhase[i] == 'stationary' and (S[t] / (k_r + S[t])) > 0.01 ):
					growthPhase[i] = 'lag'
					growthCycleTimer[i] = 0.0
				elif( growthPhase[i] == 'lag' and (rho[i]-r[i])/rho[i] < 0.5 ):
					growthPhase[i] = 'exponential'

			l[t+1] 		= lamda * ( A[t]**eta / (k_l**eta + A[t]**eta) ) * r[t]

			d_B[t+1]	= (Vmax_B * A[t])/(k_B + A[t])

			d_Bext[t+1]	= (Vmax_Bext * Aext[t])/(k_Bext + Aext[t])


		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# SAVE OUT FINAL VARIABLE VALUES: 
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		for i, inoculum in enumerate(inoculums):

			inoculum.cellCount 	= N[-1][i]
			inoculum.growthRate = r[-1][i]
			inoculum.lysisRate 	= l[-1][i]

			for solute in inoculum.periplasm.solutes:
				if(isinstance(solute, BetaLactamase)):
					solute.concentration 	= B[-1][i]
					solute.hydrolysisRate 	= d_B[-1][i]

					# Find the corresponding extracellular BetaLactamase:
					for solute_ext in media.solutes:
						if(isinstance(solute_ext, BetaLactamase) and solute == solute_ext):
							solute_ext.concentration 	= Bext[-1][i]
							solute_ext.hydrolysisRate 	= d_Bext[-1][i]
							break

					# Assume Media objects hold only one solute object per type (BetaLactamase)
					break

			for drug in inoculum.periplasm.drugs:
				if(isinstance(drug, BetaLactam)):
					drug.concentration 	= A[-1][i]
					# Assume Media objects hold only one drug object per type (BetaLactam)
					break

		#--------------------------------------------------

		for drug in media.drugs:
			if(isinstance(drug, BetaLactam)):
				drug.concentration = Aext[-1]
				# Assume Media objects hold only one drug object per type (BetaLactam)
				break

		media.nutrient.concentration = S[-1] / media.volume


		for i, strainName in enumerate([inoculum.strain.name for inoculum in inoculums]):
			t_vals 			= numpy.arange(start=0, stop=time+dt, step=dt)
			t_vals_macro 	= t_vals if self.num_runs==0 else t_vals+self.data['t'].max() 	# in terms of an overall timer, assume that each run starts immediately where the previous run left off
			run_timeseries 	= {
								'run': 	[self.num_runs]*(T+1),
								'strain': 	[strainName]*(T+1),
								't': 		t_vals_macro,
								't_run':	t_vals,
								'N': 		N[:,i],
								'r': 		r[:,i],
								'l': 		l[:,i],
								'B': 		B[:,i],
								'Bext': 	Bext[:,i],
								'd_B': 		d_B[:,i],
								'd_Bext': 	d_Bext[:,i],
								'A': 		A[:,i],
								'Aext': 	Aext[:,i],
								'S': 		S[:,i]
							  }




			strain_run_data = pandas.DataFrame.from_dict(run_timeseries)

			self.data = self.data.append(strain_run_data)
			self.data = self.data.drop_duplicates(subset=['t', 'strain'], keep='first')


		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		self.num_runs 	+= 1

		return

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# UPDATE THE DATAFRAME HOLDING THE VARIABLE TIME SERIES:
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# TODO: make this return a DataFrame instead

		output_dt = downsample_output_dt if downsample_output_dt is not None else dt
		
		timeseries = {	
						'N': 		N[::output_dt],
						'r': 		r[::output_dt],
						'l': 		l[::output_dt],
						'B': 		B[::output_dt],
						'Bext': 	Bext[::output_dt],
						'd_B': 		d_B[::output_dt],
						'd_Bext': 	d_Bext[::output_dt],
						'A': 		A[::output_dt],
						'Aext': 	Aext[::output_dt],
						'S': 		S[::output_dt]
					 }

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~s

		self.num_runs 	+= 1

		return 


	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	def figure(self):
		import matplotlib.pyplot as pyplot
		import seaborn

		fig = pyplot.figure(figsize=(16,9))
		ax1 = pyplot.subplot(221)
		ax2 = pyplot.subplot(223)
		ax3 = pyplot.subplot(222)

		seaborn.lineplot(data=self.data, ax=ax1, x="t", y="N", hue="strain")

		seaborn.lineplot(data=self.data, ax=ax2, x="t", y="r", hue="strain")
		seaborn.lineplot(data=self.data, ax=ax2, x="t", y="l", hue="strain")

		seaborn.lineplot(data=self.data, ax=ax3, x="t", y="S")

		seaborn.set_style("ticks")
		seaborn.despine()

		fig.tight_layout()

		pyplot.show()




			
			
			


		




