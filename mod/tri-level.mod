/*********************************************
 * OPL 12.10.0.0 Model
 * Author: yangk
 * Creation Date: Aug 25, 2020 at 10:45:21 AM
 *********************************************/
 
/* Tuples */
// Tuple that represents the demand attributes
tuple demandAttrTuple{
  	int i;        // Starting region
  	int j;        // Destination region
  	int t;        // Time
  	float v;      // Volume of the demand
  	int tt;       // Transfer time
  	float p;      // Price
}

// Tuple that represents the edge attributes
tuple edgeAttrTuple{
  	int i;        // Starting region
  	int j;        // Destination region
  	int t;        // Time
}

// Tuple that represents an edge
tuple Edge{
  int i;          // Starting region
  int j;          // Destination region
}

// Tuple that represents the demand edge
tuple demandEdgeTuple{
  int i;          // Starting region
  int j;          // Destination region
  int t;          // Time
}

// Tuple that represents the accumulation
tuple accTuple{
  int i;          // Region identifier
  float n;        // Number of accumulated units
}

// Tuple that represents the daily accumulation
tuple daccTuple{
  int i;          // Region identifier
  int t;          // Time
  float n;        // Number of accumulated units per day
}

/* Model Parameters */
string path = ...;                   // File path for output
int t0 = ...;                        // Start time
int T = ...;                         // End time
int tf = t0+T;                       // Final time
float beta = ...;                    // Beta value (used in objective function)
{demandAttrTuple} demandAttr = ...;  // Set of demand attributes
{edgeAttrTuple} edgeAttr = ...;      // Set of edge attributes
{accTuple} accInitTuple = ...;       // Set of initial accumulations
{daccTuple} daccAttr = ...;          // Set of daily accumulations

// Derived Sets
{Edge} edge = {<i,j>|<i,j,t> in edgeAttr};                       // Set of edges
{int} region = {i|<i,v> in accInitTuple};                        // Set of regions
float accInit[region] = [i:v|<i,v> in accInitTuple];             // Mapping of initial accumulations
float dacc[region][t0..tf-1] = [i:[t:v]|<i,t,v> in daccAttr];    // Mapping of daily accumulations


{demandEdgeTuple} demandEdge = {<i,j,t>|<i,j,t,v,tt,p> in demandAttr}; 		// Set of demand edges
float demand[demandEdge] = [<i,j,t>:v|<i,j,t,v,tt,p>in demandAttr];   		// Mapping of demand volumes
float demandTime[demandEdge] = [<i,j,t>:tt|<i,j,t,v,tt,p> in demandAttr];   // Mapping of transfer times
float price[demandEdge] = [<i,j,t>:p|<i,j,t,v,tt,p> in demandAttr];  		// Mapping of prices
int tt[edge] = [<i,j>:t|<i,j,t> in edgeAttr];                       		// Mapping of transfer times

/* Decision Variables */
dvar float+ demandFlow[edge][t0..tf-1];       // Demand flow over edges
dvar float+ rebFlow[edge][t0..tf-1];          // Rebalancing flow over edges
dvar float+ acc[region][t0..tf];              // Accumulation over time in each region
dvar float+ rho[region][t0..tf-1];            // Auxiliary variable
dvar float+ pi[region][t0..tf-1];             // Auxiliary variable
dvar float+ xi[edge][t0..tf-1];               // Auxiliary variable
dvar float+ desiredAcc[region][t0..tf];       // Desired accumulation
dvar boolean u[region][t0..tf-1];             // Boolean variable for constraints
dvar boolean v[region][t0..tf-1];             // Boolean variable for constraints
dvar boolean w[edge][t0..tf-1];               // Boolean variable for constraints
dvar float+ alpha[edge][t0..tf-1];            // Auxiliary variable
dvar float+ phi[region][t0..tf-1];            // Auxiliary variable
dvar float+ gamma[edge][t0..tf-1];            // Auxiliary variable
dvar boolean x[edge][t0..tf-1];               // Boolean variable for constraints
dvar boolean y[region][t0..tf-1];             // Boolean variable for constraints
dvar boolean z[edge][t0..tf-1];               // Boolean variable for constraints


// Objective Function

// Maximizes profit from demand flow minus costs from rebalancing and demand time.
maximize(sum(e in demandEdge) demandFlow[<e.i,e.j>][e.t]*price[e] 
       - beta * sum(e in edge,t in t0..tf-1)rebFlow[e][t]*tt[e]
       - beta * sum(e in edge,t in t0..tf-1:<e.i,e.j,t> in demandEdge)demandFlow[e][t]*demandTime[<e.i,e.j,t>]);

// Constraints

subject to
{
  // Flow conservation and bounds for each region and time.
  forall(t in t0..tf-1)
  {
    forall(i in region)
    {
      acc[i][t+1] == acc[i][t] 
              - sum(e in edge: e.i==i)(demandFlow[e][t] + rebFlow[e][t]) 
              + sum(e in demandEdge: e.j==i && e.t+demandTime[e]==t)demandFlow[<e.i,e.j>][e.t] 
              + sum(e in edge: e.j==i && t-tt[e]>=t0)rebFlow[e][t-tt[e]] 
              + dacc[i][t];
      sum(e in edge: e.i==i)(demandFlow[e][t]+ rebFlow[e][t]) <= acc[i][t];
      acc[i][t] - sum(e in edge: e.i==i)demandFlow[e][t] <= 20000 * (1-y[i][t]);
      phi[i][t] <= 20000*y[i][t];
      sum(e in edge: e.i==i && e.i!=e.j) (rebFlow[<e.j, e.i>][t] - rebFlow[<e.i, e.j>][t]) - desiredAcc[i][t] + acc[i][t]  - sum(e in edge: e.i==i)demandFlow[e][t] >=0;
      sum(e in edge: e.i==i && e.i!=e.j) (rebFlow[<e.j, e.i>][t] - rebFlow[<e.i, e.j>][t]) - desiredAcc[i][t] + acc[i][t]  - sum(e in edge: e.i==i)demandFlow[e][t]<= 20000 * u[i][t];
      acc[i][t] - sum(e in edge: e.i==i)rebFlow[e][t]  - sum(e in edge: e.i==i)demandFlow[e][t] <= 20000 * v[i][t];
      rho[i][t] <= 20000 * (1-u[i][t]);
      pi[i][t] <= 20000 * (1-v[i][t]);
      if(t == t0)
        acc[i][t] == accInit[i];
    }

    // Flow conservation and bounds for each edge and time.
    forall(e in edge)
    {
      if(<e.i,e.j,t> in demandEdge)
      {
        demandFlow[e][t] <= demand[<e.i,e.j,t>];
        -price[<e.i,e.j,t>] + alpha[e][t] + phi[e.i][t] - gamma[e][t] == 0;
        demand[<e.i,e.j,t>] - demandFlow[e][t] <= demand[<e.i,e.j,t>]*x[e][t];
      }
      else
      {
        demandFlow[e][t] == 0;      
        alpha[e][t] + phi[e.i][t] - gamma[e][t] == 0;
      }
      demandFlow[e][t] <= 20000*z[e][t];
      alpha[e][t] <= 20000*(1-x[e][t]);
      gamma[e][t] <= 20000*(1-z[e][t]);
      tt[e] - rho[e.j][t] + rho[e.i][t] + pi[e.i][t] - xi[e][t] == 0;
      rebFlow[e][t] <= 20000 * (1-w[e][t]);
      xi[e][t] <= 20000 * w[e][t];
    }
  }
}

// Main Execution

main {
  // Generates the model.
  thisOplModel.generate();
  
  // Solves the model.
  cplex.solve();
  
  // Outputs the solution to a file.
  var t = thisOplModel.t0
  var ofile = new IloOplOutputFile(thisOplModel.path);
  ofile.write("flow=[")
  for(var e in thisOplModel.edge)
    if(thisOplModel.demandFlow[e][t]>1e-3 || thisOplModel.rebFlow[e][t]>1e-3)
    {
      ofile.write("(");
      ofile.write(e.i);
      ofile.write(",");
      ofile.write(e.j);
      ofile.write(",");
      ofile.write(thisOplModel.demandFlow[e][t]);
      ofile.write(",");
      ofile.write(thisOplModel.rebFlow[e][t]);
      ofile.write(")");
    }
  ofile.writeln("];")
  ofile.write('desiredAcc=[')
  for(var i in thisOplModel.region)
  {
    ofile.write("(");
    ofile.write(i);
    ofile.write(",");
    ofile.write(thisOplModel.desiredAcc[i][t]);
    ofile.write(")");
  }
  ofile.writeln("];")
  ofile.close();
}