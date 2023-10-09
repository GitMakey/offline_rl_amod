/*********************************************
 * OPL 12.10.0.0 Model
 * Author: yangk
 * Creation Date: Aug 25, 2020 at 10:45:21 AM
 * Purpose: Optimizing flow of entities between regions over time, maximizing revenue from demand fulfillment while balancing costs.
 *********************************************/

// --- TUPLE DEFINITIONS ---

// Describes demand attributes for two regions at a specific time.
tuple demandAttrTuple{
  	int i;      // Start region
  	int j;      // End region
  	int t;      // Specific time of demand
  	float v;    // Volume of demand
  	float tt;   // Travel time between i and j
  	float p;    // Price or revenue potential for this demand
}

// Describes a connection or edge between two regions at a specific time.
tuple edgeAttrTuple{
  	int i;
  	int j;
  	int t;
}

// Simplified connection between two regions, without time specification.
tuple Edge{
  int i;
  int j;
}

// Represents demand between two regions at a specific time.
tuple demandEdgeTuple{
  int i;
  int j;
  int t;
}

// Describes the current accumulation or stock of entities in a region.
tuple accTuple{
  int i;
  float n;
}

// Describes the change in accumulation in a region over time.
tuple daccTuple{
  int i;
  int t;
  float n;
}

// --- INPUT DATA ---

string path = ...;            // Path for the output file
int t0 = ...;                 // Starting time
int T = ...;                  // Total time period
int tf = t0+T;                // Calculated end time
float beta = ...;             // Cost factor
{demandAttrTuple} demandAttr = ...;
{edgeAttrTuple} edgeAttr = ...;
{accTuple} accInitTuple = ...;
{daccTuple} daccAttr = ...;

// --- DERIVED DATA ---

// Extracting the simple edge data from the edge attributes
{Edge} edge = {<i,j>|<i,j,t> in edgeAttr};
{int} region = {i|<i,v> in accInitTuple}; 
float accInit[region] = [i:v|<i,v> in accInitTuple];
float dacc[region][t0..tf-1] = [i:[t:v]|<i,t,v> in daccAttr];
{demandEdgeTuple} demandEdge = {<i,j,t>|<i,j,t,v,tt,p> in demandAttr};

// Mapping demand attributes to arrays for easy access.
float demand[demandEdge] = [<i,j,t>:v|<i,j,t,v,tt,p> in demandAttr];
float price[demandEdge] = [<i,j,t>:p|<i,j,t,v,tt,p> in demandAttr];
float demandTime[demandEdge] = [<i,j,t>:tt|<i,j,t,v,tt,p> in demandAttr];
int tt[edge] = [<i,j>:t|<i,j,t> in edgeAttr];

// --- DECISION VARIABLES ---

dvar float+ demandFlow[edge][t0..tf-1];  // Flow for demand between regions over time
dvar float+ rebFlow[edge][t0..tf-1];     // Rebalancing flow between regions over time
dvar float+ acc[region][t0..tf];         // Accumulation in regions over time

// --- OBJECTIVE FUNCTION ---

// Maximize revenue considering price, rebalancing cost, and demand fulfillment cost
maximize(sum(e in demandEdge) demandFlow[<e.i,e.j>][e.t]*price[e] 
         - beta * sum(e in edge,t in t0..tf-1)rebFlow[e][t]*tt[e] 
         - beta * sum(e in edge,t in t0..tf-1:<e.i,e.j,t> in demandEdge)demandFlow[e][t]*demandTime[<e.i,e.j,t>]);

// --- CONSTRAINTS ---

subject to
{
  forall(t in t0..tf-1)
  {
    // Ensure flow conservation and update accumulation correctly
    forall(i in region)
    {  
    	acc[i][t+1] == acc[i][t] - sum(e in edge: e.i==i)(demandFlow[e][t] + rebFlow[e][t]) 
      			+ sum(e in demandEdge: e.j==i && e.t+demandTime[e]==t)demandFlow[<e.i,e.j>][e.t] + sum(e in edge: e.j==i && t-tt[e]>=t0)rebFlow[e][t-tt[e]] + dacc[i][t];
		sum(e in edge: e.i==i)(demandFlow[e][t]+ rebFlow[e][t]) <= acc[i][t];
      	if(t == t0)
      		acc[i][t] == accInit[i];
 	}  	    
    // Ensure demand is met and flows are correctly initialized
    forall(e in edge)
      if(<e.i,e.j,t> in demandEdge)
      		demandFlow[e][t] <= demand[<e.i,e.j,t>];
      else
      		demandFlow[e][t] == 0;      
  }
  
}

// --- MAIN EXECUTION ---

main {
  thisOplModel.generate();
  cplex.solve();

  // Output results to file
  var t = thisOplModel.t0;
  var ofile = new IloOplOutputFile(thisOplModel.path);
  ofile.write("flow=[");
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
  ofile.close();
}