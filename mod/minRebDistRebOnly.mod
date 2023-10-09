/**
 * Description:
 * The script aims to minimize the rebalancing time of vehicles across regions while satisfying certain conditions.
 * The model assumes each region has a certain initial number of vehicles and desired number of vehicles.
 * The model determines the best rebalancing flow across edges to minimize the rebalancing time.
 * 
 * Key Components:
 * - Edges representing possible routes for rebalancing.
 * - Time taken on each edge.
 * - Initial and desired number of vehicles in each region.
 * - Decision variables for the rebalancing flow.
 * - Objective function to minimize the total rebalancing time.
 * - Constraints to ensure proper vehicle rebalancing.
 */

// Set CPLEX to use 4 threads.
execute
{
  cplex.threads=4;
}

// Tuple definition for an edge.
tuple Edge{
  int i;      // Starting node of the edge.
  int j;      // Ending node of the edge.
}

// Tuple definition for edge attributes.
tuple edgeAttrTuple{
    int i;   // Starting node of the edge.
    int j;   // Ending node of the edge.
    int t;   // Time taken on the edge.
}

// Tuple definition for account initialization.
tuple accTuple{
  int i;     // Region index.
  float n;   // Number of vehicles.
}

// Path for saving the solution.
string path = ...;

// Input data for edge attributes, initial account, and desired vehicles.
{edgeAttrTuple} edgeAttr = ...;
{accTuple} accInitTuple = ...;
{accTuple} accRLTuple = ...;

// Construct edge set from edgeAttr data.
{Edge} edge = {<i,j>|<i,j,t> in edgeAttr};

// Construct region set from accInitTuple data.
{int} region = {i|<i,v> in accInitTuple};

// Create mappings for time taken on each edge and the desired and initial number of vehicles in each region.
float time[edge] = [<i,j>:t|<i,j,t> in edgeAttr];
float desiredVehicles[region] = [i:v|<i,v> in accRLTuple];
float vehicles[region] = [i:v|<i,v> in accInitTuple];

// Decision variables for the rebalancing flow across edges.
dvar int+ demandFlow[edge];
dvar int+ rebFlow[edge];

// Objective function to minimize the total rebalancing time.
minimize(sum(e in edge) (rebFlow[e]*time[e]));

// Constraints to ensure proper vehicle rebalancing.
subject to
{
  forall(i in region)
    {
    // Ensure the net rebalancing flow satisfies the desired number of vehicles in the region.
    sum(e in edge: e.i==i && e.i!=e.j) (rebFlow[<e.j, e.i>] - rebFlow[<e.i, e.j>]) >= desiredVehicles[i] - vehicles[i];
    // Ensure the rebalancing flow out of a region does not exceed the initial number of vehicles in the region.
    sum(e in edge: e.i==i && e.i!=e.j) rebFlow[<e.i, e.j>] <= vehicles[i];
    }
}

// Main procedure to solve the model and write the solution to a file.
main {
  thisOplModel.generate();
  cplex.solve();

  var ofile = new IloOplOutputFile(thisOplModel.path);
  ofile.write("flow=[");

  // Write the rebalancing flow for each edge to the file.
  for(var e in thisOplModel.edge)
       {
         ofile.write("(");
         ofile.write(e.i);
         ofile.write(",");
         ofile.write(e.j);
         ofile.write(",");
         ofile.write(thisOplModel.rebFlow[e]);
         ofile.write(")");
       }

  ofile.writeln("];")

  // Write the objective function value to the file.
  var obj = cplex.getObjValue();
	ofile.writeln("obj="+obj+";");

  ofile.close();
}
