package wsc.ecj.nsga2;

import static java.lang.Math.abs;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;

import ec.BreedingPipeline;
import ec.EvolutionState;
import ec.Individual;
import ec.multiobjective.MultiObjectiveFitness;
import ec.util.Parameter;
import wsc.data.pool.InitialWSCPool;
import wsc.data.pool.Service;
import wsc.graph.ServiceGraph;
import wsc.problem.WSCInitializer;

public class WSCSmartLocalSearchPipeline extends BreedingPipeline {

	private static final long serialVersionUID = 1L;

	@Override
	public Parameter defaultBase() {
		return new Parameter("wscmutationpipeline");
	}

	@Override
	public int numSources() {
		return 1;
	}

	@Override
	public int produce(int min, int max, int start, int subpopulation, Individual[] inds, EvolutionState state,
			int thread) {

		int n = sources[0].produce(min, max, start, subpopulation, inds, state, thread);

		if (!(sources[0] instanceof BreedingPipeline)) {
			inds[start] = (Individual) (inds[start].clone());
		}

		if (!(inds[start] instanceof SequenceVectorIndividual))
			// uh oh, wrong kind of individual
			state.output
					.fatal("WSCMutationPipeline didn't get a SequenceVectorIndividual. The offending individual is: "
							+ inds[start]);

		WSCInitializer init = (WSCInitializer) state.initializer;
		SequenceVectorIndividual bestNeighbour = (SequenceVectorIndividual) inds[start].clone();

		double bestScore;
		if (WSCInitializer.tchebycheff)
			bestScore = init.calculateTchebycheffScore(bestNeighbour, start);
		else
			bestScore = init.calculateScore(bestNeighbour, start);

		SequenceVectorIndividual neighbour;

		int localSearchTries = 0;
		for (int indexA = 0; indexA < bestNeighbour.genome.length; indexA++) {
			for (int indexB = 0; indexB < bestNeighbour.genome.length; indexB++) {
				if (localSearchTries >= WSCInitializer.numLocalSearchTries)
					break;

				neighbour = (SequenceVectorIndividual) inds[start].clone();

				Set<Service> usedServices = Sets.newHashSet(neighbour.genome);

				if (swapIsAllowed(neighbour.genome, indexA, indexB, usedServices, start, init)) {
					swapServices(neighbour.genome, indexA, indexB);
					localSearchTries++;

					// Calculate fitness of neighbor
					neighbour.calculateSequenceFitness(neighbour, init, state);

					double score;
					if (WSCInitializer.tchebycheff)
						score = init.calculateTchebycheffScore(neighbour, start);
					else
						score = init.calculateScore(neighbour, start);

					// If the neighbour has a better fitness score than the current best, set
					// current best to be neighbour
					if (score < bestScore) {
						bestScore = score;
						bestNeighbour = neighbour;
					}
				}
			}
		}

		inds[start] = bestNeighbour;
		inds[start].evaluated = false;

		return n;
	}

	private boolean swapIsAllowed(Service[] genome, int indexA, int indexB, Set<Service> usedServices, int problemIndex,
			WSCInitializer init) {
		boolean swap = false;
		// Do not swap services if both of them are unused. Only swap if one is used and
		// the other is not.
		if (usedServices.contains(genome[indexA]) && !usedServices.contains(genome[indexB])) {
//			Service serviceA = genome[indexA];
//			Service serviceB = genome[indexB];
//			boolean better = qualityIsBetter(serviceA, serviceB, problemIndex, init);
//			// Only perform swap if used service does not completely dominate the unused
//			// one.
//			if (!better)
				swap = true;

		} else if (usedServices.contains(genome[indexB]) && !usedServices.contains(genome[indexA])) {
			// Do not swap services if the QoS of the unused service is lower than that of
			// the used one
			// Service serviceA = genome[indexA];
			// Service serviceB = genome[indexB];
			// boolean better = qualityIsBetter(serviceB, serviceA, problemIndex, init);
			// // Only perform swap if used service does not completely dominate the unused
			// // one.
			// if (!better)
				swap = true;
		}

		return swap;
	}

	private void swapServices(Service[] genome, int indexA, int indexB) {
		Service temp = genome[indexA];
		genome[indexA] = genome[indexB];
		genome[indexB] = temp;
	}

//	private boolean qualityIsBetter(Service serviceA, Service serviceB, int problemIndex, WSCInitializer init) {
//		// Calculate Tchebycheff quality for individual services
//		double scoreA = calculateIndividualTchebycheffScore(serviceA, problemIndex, init);
//		double scoreB = calculateIndividualTchebycheffScore(serviceB, problemIndex, init);
//
//		return scoreA < scoreB;
//	}

//	public double calculateIndividualTchebycheffScore(Service serv, int problemIndex, WSCInitializer init) {
//		double[] problemWeights = init.weights[problemIndex];
//		double max_fun = -1 * Double.MAX_VALUE;
//
//		// Calculate objectives for service
//		double[] qos = serv.getQos();
//
//		double a = normaliseIndividualAvailability(qos[WSCInitializer.AVAILABILITY], init);
//		double r = normaliseIndividualReliability(qos[WSCInitializer.RELIABILITY], init);
//		double t = normaliseIndividualTime(qos[WSCInitializer.TIME], init);
//		double c = normaliseIndividualCost(qos[WSCInitializer.COST], init);
//
//		double[] objectives = new double[2];
//		objectives[0] = t + c;
//		objectives[1] = a + r;
//
//		for (int i = 0; i < problemWeights.length; i++) {
//			double diff = abs(objectives[i] - init.idealPoint[i]);
//			double feval;
//			if (problemWeights[i] == 0)
//				feval = 0.00001 * diff;
//			else
//				feval = problemWeights[i] * diff;
//			if (feval > max_fun)
//				max_fun = feval;
//		}
//		return max_fun;
//	}

//	private double normaliseIndividualAvailability(double availability, WSCInitializer init) {
//		if (init.maxIndividualAvailability - init.minIndividualAvailability == 0.0)
//			return 1.0;
//		else
//			// return (availability - init.minAvailability)/(init.maxAvailability -
//			// init.minAvailability);
//			return (init.maxIndividualAvailability - availability)
//					/ (init.maxIndividualAvailability - init.minIndividualAvailability);
//	}
//
//	private double normaliseIndividualReliability(double reliability, WSCInitializer init) {
//		if (init.maxIndividualReliability - init.minIndividualReliability == 0.0)
//			return 1.0;
//		else
//			// return (reliability - init.minReliability)/(init.maxReliability -
//			// init.minReliability);
//			return (init.maxIndividualReliability - reliability)
//					/ (init.maxIndividualReliability - init.minIndividualReliability);
//	}
//
//	private double normaliseIndividualTime(double time, WSCInitializer init) {
//		if (init.maxIndividualTime - init.minIndividualTime == 0.0)
//			return 1.0;
//		else
//			// return (init.maxTime - time)/(init.maxTime - init.minTime);
//			return (time - init.minIndividualTime) / (init.maxIndividualTime - init.minIndividualTime);
//	}
//
//	private double normaliseIndividualCost(double cost, WSCInitializer init) {
//		if (init.maxIndividualCost - init.minIndividualCost == 0.0)
//			return 1.0;
//		else
//			// return (init.maxCost - cost)/(init.maxCost - init.minCost);
//			return (cost - init.minIndividualCost) / (init.maxIndividualCost - init.minIndividualCost);
//	}
}
