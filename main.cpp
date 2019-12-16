#include <cmath>
#include <vector>
#include <utility>
#include <iostream>
#include <random>
#include <iomanip>
#include <algorithm>

constexpr size_t PopulationSize = 4;

constexpr double MutationPossibility = 0.25;
constexpr double MutationShift = 0.1;

constexpr std::pair<double, double> XInterval = std::make_pair(-2, 2);
constexpr std::pair<double, double> YInterval = std::make_pair(-2, 2);

double FitnessFunction(std::pair<double, double> point) {
    return -log(1 + pow(point.first, 2) + pow(point.second, 2)) + 3;
}

typedef std::pair<double, double> Person;
typedef std::vector<Person> Population;

void PrintHeader();

Population SourcePopulation();

void PrintTable(size_t generation, const Population& population);
double AverageFit(const Population& population);
double MaximumFit(const Population& population);

Population Selection(const Population& population);
Population Crossover(const Population& sorted);
Population Mutation(Population population);

double MutateGen(double gen);

double Random(double a, double b);

int main() {
    PrintHeader();

    Population population = SourcePopulation();
    PrintTable(0, population);

    for (size_t generation = 1; generation <= 100; ++generation) {
        population = Mutation (Crossover (Selection(population)));

        PrintTable(generation, population);
    }

    return 0;
}

double Random(double a, double b) {
    static std::random_device rd;
    std::uniform_real_distribution<double> distribution(a, b);

    return distribution(rd);
}

void PrintHeader() {
    std::cout << "___________________________________________________________________" << std::endl
              << "|Generation|     X    |     Y    |    Fit   |    Max   |    Ave   |" << std::endl;
}

Population SourcePopulation() {
    Population population;
    for (size_t i = 0; i < PopulationSize; ++i) {
        population.emplace_back(std::make_pair(Random(XInterval.first, XInterval.second),
                                               Random(YInterval.first, YInterval.second)
                                               )
                                );
    }

    return population;
}

void PrintTable(size_t generation, const Population& population) {
    bool isJointPrint = false;
    for (const auto& person : population) {
        if (!isJointPrint) {
            std::cout << "|" << std::setw(10) << generation;
        } else {
            std::cout << "|" << std::setw(10) << "";
        }

        std::cout << "|" << std::setw(10) << std::setprecision(5) << person.first
                  << "|" << std::setw(10) << std::setprecision(5) << person.second
                  << "|" << std::setw(10) << std::setprecision(5) << FitnessFunction(person);

        if (!isJointPrint) {
            isJointPrint = true;

            std::cout << "|" << std::setw(10) << std::setprecision(5) << MaximumFit(population)
                      << "|" << std::setw(10) << std::setprecision(5) << AverageFit(population)
                      << "|" << std::endl;
        } else {
            std::cout << "|" << std::setw(10) << ""
                      << "|" << std::setw(10) << "" << "|"
                      << std::endl;
        }
    }
    std::cout << std::endl;
}

double SumFit(const Population& population) {
    return std::accumulate(population.begin(), population.end(), 0.,
                           [](double& bucket, auto person) {
                               return bucket += FitnessFunction(person);
                           });
}

double AverageFit(const Population& population) {
    return SumFit(population) / population.size();
}

double MaximumFit(const Population& population) {
    auto found = *std::max_element(population.begin(),
                                   population.end(),
                                   [](auto first, auto second) {
                                       return FitnessFunction(first) < FitnessFunction(second);
                                   }
                                   );

    return FitnessFunction(found);
}



std::vector<double> MakeRoulette(const Population& population) {
    std::vector<double> roulette;
    if (population.empty()) {
        return roulette;
    }

    double bound = 0;
    for (const auto& person : population) {
        bound += (FitnessFunction(person)) /
                  SumFit(population);
        roulette.push_back(bound);
    }

    return roulette;
}

size_t Shot(const std::vector<double>& roulette) {
    auto shot = Random(0., 1.);

    auto iter = std::lower_bound(roulette.begin(), roulette.end(), shot);

    return std::distance(roulette.begin(), iter);
}

Population Selection(const Population& population) {
    auto roulette = MakeRoulette(population);

    Population sorted;
    for (const auto& person : population) {
        sorted.emplace_back(population[Shot(roulette)]);
    }

    return sorted;
}

Population Crossover(const Person& first, const Person& second) {
    Population children;
    children.emplace_back(std::make_pair(first.first, second.second));
    children.emplace_back(std::make_pair(second.first, first.second));

    return children;
}

Population Crossover(const Population& sorted) {
    Population population;

    Population children = Crossover(sorted[0], sorted[1]);
    population.insert(population.end(), children.begin(), children.end());

    children = Crossover(sorted[0], sorted[2]);
    population.insert(population.end(), children.begin(), children.end());

    return population;
}

Population Mutation(Population population) {
    for (auto& person : population) {
        if (Random(0., 1.) < MutationPossibility) {
            person.first = MutateGen(person.first);
        }

        if (Random(0., 1.) < MutationPossibility) {
            person.second = MutateGen(person.second);
        }
    }

    return population;
}

double MutateGen(double gen) {
    if (Random(0., 1.) < 0.5) {
        if (gen - MutationShift > XInterval.first) {
            gen -= MutationShift;
        }/* else {
            gen += MutationShift;
        }*/
    } else {
        if (gen + MutationShift < XInterval.second) {
            gen += MutationShift;
        }/* else {
            gen -= MutationShift;
        }*/
    }

    return gen;
}
