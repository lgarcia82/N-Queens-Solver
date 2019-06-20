/*
	Title: N Queens Puzzle
	Description: This application implements a genetic algorithm to find and verify all solutions
	for the N queen puzzle for all values of N, such as 0 < N < 10. All methods of genetic algorthims 
	are implemented and utilized.
	This applicaton implements and utilizes a fitness function, selection operator, crossover operator, and mutation operator
	to achieve all desired results
	Source file: Main.cpp
	Data file: NA
	Authors: Tri Vo & Lee Garcia
	Date: Spring 2018
*/
#include<iostream>
#include<cmath>
#include<cstdlib>
#include<algorithm>
#include<sstream>
#include<cstring>
#include<map>
#include<list>
#include<queue>
#include<ctime>

using namespace std;
//population node
typedef struct
{
	string arrangement;
	int cost;
	int iteration;
} individual;

//vector data structure of node type
typedef vector<individual*> population_type;
//map data structure to map nodes which meet fitness goal and are classified as a solution
typedef map<string, int> SolutionsFound;

population_type population; //initial vector object
int chessBoardSize;          //user defined attribute which defines board size and # of queens
int initialPopulationCount = 1000;  //initial population size of each generation before breeding

/*
	fitness function
	@param: arrangement member of individual node
	captures the fitness score of each individual. Works in a decrementing manner. 
	equation for ideal fitness => (N * (N - 1))/2
*/
int fitnessValue(string arrangement)
{
	int fitness = (chessBoardSize*(chessBoardSize - 1)) / 2;          //establish ideal fitness score according to N's value
																	  //ideal fitness for N = 8 is 28
	for (int i = 0; i<chessBoardSize; i++)
		for (int j = i + 1; j<chessBoardSize; j++)    //removing pairs that lie on the same row and on the same diagonal
			if (/*(arrangement[i] == arrangement[j]) || */(i - arrangement[i] == j - arrangement[j]) || (i + arrangement[i] == j + arrangement[j]))
				fitness--;      //test each queens position in the string arrangement for (+- 1) condition. 
								//every flag reduces fitness number by 1
	return fitness;
}

/*
	function to create individual node
	@param: none
	creates single nodes of individual type
*/
individual* createNode()
{
	individual *newNode = new individual;
	return newNode;
}

/*
	function to create initial population before breeding process
	@param: map data structure of found solutions
	generates initial population and every sequential generation that follows
*/
void generatePopulation(SolutionsFound set)
{
	string sampleArrangement = "";		//arrangement declaration
	individual *temp;					// temp individual for population generation

	//create a string arrangement of (0 - N), only distinct values exist in arrangement
	for (int i = 0; i < chessBoardSize; i++)
	{
		ostringstream ostr;
		ostr << i;
		sampleArrangement += ostr.str();
	}

	//generate population by shuffling sample arrangement and storing in population vector
	int fitnessGoal = (chessBoardSize*(chessBoardSize - 1)) / 2;
	for (int i = 0; i<initialPopulationCount; i++)
	{
		random_shuffle(sampleArrangement.begin(), sampleArrangement.end());
		temp = createNode();
		temp->arrangement = sampleArrangement;
		temp->cost = fitnessValue(sampleArrangement); // pass arrangement into fitness function
		population.push_back(temp); // store arrangement into population
		
	}
	
	// take each generated individual and evaluate cost 
	for (int i = 0; i < population.size(); i++)  //eliminate unfit chromosomes
	{
		if (population[i]->cost < fitnessGoal)
		{
			population.erase(population.begin() + i); //if cost of individual is NOT ideal, delete from the population
		}
	}

	//create map iterator
	SolutionsFound::iterator pos;

	//search population for node arrangements that exist as solutions already found
	for (pos = set.begin(); pos != set.end(); ++pos)   //remove found solution from the population
	{
		for (int i = 0; i < population.size(); i++)
		{
			if (pos->first == population[i]->arrangement)   //if arrangement exist in found solution map, delete from population
			{
				population.erase(population.begin() + i);
			}
		}
	}
}

/*
	crossover operator, single point operation
	@param: individual node x, individual node y - parents 
*/
individual* reproduce(individual *x, individual *y)
{
	individual *child = createNode(); // create node to hold child
	int n = chessBoardSize;
	int c = rand() % n;    //pick random point for crossover from 0 - N
	child->arrangement = (x->arrangement).substr(0, c) + (y->arrangement).substr(c, n - c + 1); // combine parent x and parent y at crossover point indicated
	child->cost = fitnessValue(child->arrangement); // evaluate fitness of child
	return child;
}
/*
	 mutation operator, single mutation
	 @param: child node z
*/

individual* mutate(individual *child)
{
	int point1 = rand() % (chessBoardSize); //pick 1st random point in string index of child arrangement
	int point2 = rand() % (chessBoardSize); //pick 2nd random point in string index of child arrangement
	char temp;
	temp = child->arrangement[point1];                           //pull integer values at both points and swap
	child->arrangement[point1] = child->arrangement[point2];
	child->arrangement[point2] = temp;
	return child;                             // return mutated child
}

/*
	selection operator
	@param: none
	chooses 2 nodes as parents for breeding, picks random node from population at any index from 0 - 99
*/
int randomSelection()
{
	int randomPos = rand() % 100;
	return randomPos;
}

/*
	fitness operator
	@param: individual node from population
	runs fitness function to test child for fitness
*/
bool isFit(individual *test)
{
	if (fitnessValue(test->arrangement) == ((chessBoardSize*(chessBoardSize - 1)) / 2)) // is child ideal fitness
		return true;
	return false;
}

/*
	comparison operator used to sort population by fitness
	@param: any 2 neighboring nodes in the population
*/
bool comp(individual *a, individual*b)
{
	return(a->cost > b->cost);
}

/*
	main algorithm function
	@param: none
	primary algorithm where the genetic process occurs. In this function parent nodes 
	are chosen, breeded, and children are produced and mutated
*/
individual* GA()
{
	int randomNum1, randomNum2;					//2 random numbers for parent selection
	individual *individualX, *individualY, *child;		//parent nodes & children node
	child = NULL;
	bool found = 0;          //boolean to terminate function between solutions
	int count = 0;			// iteration counter
	while (!found)
	{
		population_type new_population;                //new population vector of individual type to hold children 
		for (unsigned int i = 0; i<population.size(); i++)						
		{
			count++;																		
			sort(population.begin(), population.end(), comp);				//sort initial population by fitness score

			randomNum1 = randomSelection();								//pick parent X
			individualX = population[randomNum1];

			randomNum2 = randomSelection();								//pick parent Y
			individualY = population[randomNum2];

			child = reproduce(individualX, individualY);				//create child in crossover

			if (rand() % 100 == 0)     //random probability for mutation of 1%
				child = mutate(child);

			if (isFit(child))         //test child for fitness
			{
				found = 1;
				child->iteration = count;
				return child;				// if child is of ideal fitness end while, return to main
			}
			new_population.push_back(child);  //else add child to new population vector
		}
		population = new_population;       // if no solution is found in 1st generation breeding, repeat process to second generation
	}
	return child;
}

/*
	initialize board size
	@param: user input
*/
void initialize(int input)
{
	srand(time(0));     //to ensure perfect randomness
	chessBoardSize = input;
}

int main()
{
	int input;                          //user input variable for board size
	long totalIterations = 0;			//variable for total iterations
	float avgIterations = 0;			//variable for average iterations per solution

	cout << "Please enter the number of queens: ";
	cin >> input;
	initialize(input);
	int maxSolutions[10] = { 1,0,0,2,10,4,40,92,352,724 };  //array with max solution based on N size. terminates program once solution max is reached
	clock_t start_time, end_time;           //to keep a track of the time spent in computing
	SolutionsFound set;						//data structure of type map to store found solutions
	int numFound = 0;						//solution counter
	int x = 1;								// map index variable
	start_time = clock();					//get start time from clock function
	if (maxSolutions[input - 1] == 0)		// no solutions
		cout << "No solution" << endl;
	else
	{
		cout << "*Returns the column number corresponding to the row at the index*" << endl << endl; 

		while (numFound != maxSolutions[input - 1])    //run until all solutions are found for N case
		{
			generatePopulation(set);					//initialize population
			individual *solution = GA();				// run genetic alogorithm and return solution when found
			if (!set[solution->arrangement])
			{
				set[solution->arrangement] = x;						//set key in map according to index
				cout << "  Solution #" << (++numFound) << ": [" << solution->arrangement << "] ---> " ; //display results every time a solution is found
				cout << "# of iterations: " << solution->iteration << endl;
				totalIterations += solution->iteration;
				x++;											//increment map value 
			}
		}
	}
	end_time = clock();      // end time for execution 
	avgIterations = totalIterations / static_cast<float>(maxSolutions[input - 1]);   //calculate average iterations per solution found
	
	// display required results
	cout << "\n\n	Execution Time: " << 1000 * ((double)(end_time - start_time) / CLOCKS_PER_SEC) << " milliseconds" << "\n";
	cout << "   Avg # of iterations: " << avgIterations << "\n\n";
	system("PAUSE");
	return 0;
}// end program Main