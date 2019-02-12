(* ::Package:: *)

(* :Package Version: 3.0.03.20  *)

(* :Copyright: Copyright 2015, Selwyn Hollis *)

(* :Name: RDL`ReactionDiffusionLab` *)
 
(* :Author: Selwyn Hollis *)

(*  New in version 1.5:
	The alternating-direction implicit method (ADI) replaces Crank-Nicolson.

	New in version 2.0:
	Support for non-square rectangular domains and periodic boundary conditions.
	Significant speedup by use of SparseArrays and LinearSolveFunctions.

	New in version 3.0:
	Ditched the surface plots; now you get only density plots.
	Dynamic output. Simplified arguments and options. Lots of general improvements.
   
 *)

If[MemberQ[Names["Global`*"], #], Remove[#]]&/@
	{"RDDensityPlots","Domain","Bounds","BoundaryConditions"};

BeginPackage["RDL`ReactionDiffusionLab`"]
 
Unprotect["RDL`ReactionDiffusionLab`*" ];
ClearAll[ "RDL`ReactionDiffusionLab`*" ];
 
RDDensityPlots::usage =
"RDDensityPlots[{f[u,v], g[u,v]}, {u,v}, {u0,v0}, {d1,d2}, {t0,t1,dt}, options]
or
RDDensityPlots[{f[u,v,w], g[u,v,w]}, {u,v,w}, {u0,v0,w0}, {d1,d2,d3}, {t0,t1,dt}, options]";

Options[RDDensityPlots]= 
	Union[{ NPlots->Automatic, Domain->{{0,Automatic},{0,1}}, Bounds->{0,1000}, Spacings->0, 
	BoundaryConditions->"Neumann", DefaultAnimationSpeedMultiplier->1.5}, Options[ListDensityPlot] ];
	
SetOptions[RDDensityPlots, PerformanceGoal->"Speed", Frame->False, Mesh->False, FrameTicks->None,
	ColorFunction->"BlueGreenYellow", PlotRange->All, AspectRatio->Automatic, ImageSize->Medium]; 

(*  Messages  *)

	RDDensityPlots::outbnds="Bounds violated. Try decreasing the time step, decreasing the size of the initial data, or widening bounds with the Bounds option.";
	RDDensityPlots::pr="Invalid PlotRange";
	RDDensityPlots::bc="Invalid or unsupported BoundaryConditions";	
	RDDensityPlots::nnfns="Can't use nonnumeric functions";
	RDDensityPlots::nndata="Can't use nonnumeric initial data";

Begin["`Private`"]
	 
(* Routines used by RDstep__ functions (defined below) *)

neumannpad= Compile[{{vals,_Real,2}},
	With[{q=(Append[#,Last[#]]&)/@(Prepend[#,First[#]]&/@vals)},
		Append[Prepend[q, First[q]], Last[q]]
		] ];

neumannXhalfstep= Compile[{{vals,_Real,2},{hh,_Real}},
	With[{kernel={{0,0,0},{1,-2,1},{0,0,0}}},
		vals + 0.5*hh*ListCorrelate[kernel, neumannpad[vals] ]
		] ];
		
neumannYhalfstep= Compile[{{vals,_Real,2},{hh,_Real}},
	With[{kernel={{0,1,0},{0,-2,0},{0,1,0}}},
		vals + 0.5*hh*ListCorrelate[kernel, neumannpad[vals] ]
		] ];
							
dirichletpad= Compile[{{vals,_Real,2}},
	Module[{q, n=Last[Dimensions[vals]], o}, o=Table[0,{n+2}];
		q=(PadLeft[#,n+2]&)/@(PadRight[#,n+1]&/@vals);
		Append[Prepend[q,o],o]
		] ];

dirichletXhalfstep= Compile[{{vals,_Real,2},{hh,_Real}},
	With[{kernel={{0,0,0},{1,-2,1},{0,0,0}}},
		vals + 0.5*hh*dirichletpad[ ListCorrelate[kernel,vals] ]
		] ];
		
dirichletYhalfstep= Compile[{{vals,_Real,2},{hh,_Real}},
	With[{kernel={{0,1,0},{0,-2,0},{0,1,0}}},
		vals + 0.5*hh*dirichletpad[ ListCorrelate[kernel,vals] ]
		] ];

setedges= Compile[{{vals,_Real,2},{bd,_Real,2},{a010,_Integer,2}}, vals*a010 + bd];

periodicize= Compile[{{vals,_Real,2}},
        With[{q=(ReplacePart[#,(First[#]+Last[#])/3 + (#[[2]]+#[[-2]])/6, {{1},{-1}} ]&)/@vals},
		ReplacePart[q, (First[q]+Last[q])/3 + (q[[2]]+q[[-2]])/6, {{1},{-1}} ]
		] ];
       
periodicXhalfstep= Compile[{{vals,_Real,2},{hh,_Real}},
	vals + 0.5*hh*(ListCorrelate[{ 1, -2, 1 }, #, {2,2} ]& /@ vals) ];
		
periodicYhalfstep= Compile[{{vals,_Real,2},{hh,_Real}},
	vals + 0.5*hh*ListCorrelate[ {{0,1,0},{0,-2,0},{0,1,0}}, vals, {2,2} ] ];

ADIXu[rhs_]:= Partition[ uLsolvex[ Flatten[rhs] ], n];

ADIYu[rhs_]:= Transpose@Partition[ uLsolvey[ Flatten[Transpose@rhs] ], m];

ADIXv[rhs_]:= Partition[ vLsolvex[ Flatten[rhs] ], n];

ADIYv[rhs_]:= Transpose@Partition[ vLsolvey[ Flatten[Transpose@rhs] ], m];

ADIXw[rhs_]:=Partition[ wLsolvex[ Flatten[rhs] ], n];

ADIYw[rhs_]:= Transpose@Partition[ wLsolvey[ Flatten[Transpose@rhs] ], m];

(** Definitions of RDstep__ functions. 
    There's one for each combination of boundary conditions. **)

RDstepNN[{u0_, v0_}, {fc_,gc_}, {d1_,d2_}, {hhx_,hhy_}]:= 
	Module[{u1,v1},
		u1= ADIXu@ (neumannYhalfstep[u0, d1*hhy ] + fc[u0,v0] );
		v1= ADIXv@ (neumannYhalfstep[v0, d2*hhy ] + gc[u0,v0] );
		{   ADIYu@ (neumannXhalfstep[u1, d1*hhx ] + fc[u1,v1] ),
		    ADIYv@ (neumannXhalfstep[v1, d2*hhx ] + gc[u1,v1] )}
		 ];
		 
RDstepPeriodic[{u0_, v0_}, {fc_,gc_}, {d1_,d2_}, {hhx_,hhy_}]:= 
	Module[{u1,v1},
		u1= ADIYu@(periodicXhalfstep[u0, d1*hhx ] + fc[u0,v0] );
		v1= ADIXv@(periodicYhalfstep[v0, d2*hhy ] + gc[u0,v0] );
	 	{  ADIXu@(periodicYhalfstep[u1, d1*hhy ] + fc[u1,v1] ),
		   ADIYv@ (periodicXhalfstep[v1, d2*hhx ] + gc[u1,v1] ) }
		 ];

RDstepDD[{u0_, v0_}, {fc_,gc_}, {d1_,d2_}, {hhx_,hhy_}]:= 
	Module[{u1,v1},
		u1= ADIXu@ (dirichletYhalfstep[u0, d1*hhy ] + fc[u0,v0] );
		v1= ADIXv@ (dirichletYhalfstep[v0, d2*hhy ] + gc[u0,v0] );
		{setedges[ ADIYu@ (dirichletXhalfstep[u1, d1*hhx ] + fc[u1,v1]), First[bdata], a010],
		 setedges[ ADIYv@ (dirichletXhalfstep[v1, d2*hhx ] + gc[u1,v1]), Last[bdata],  a010]}];
		 
RDstepND[{u0_, v0_}, {fc_,gc_}, {d1_,d2_}, {hhx_,hhy_}]:= 
	Module[{u1,v1},
		u1= ADIXu@ (  neumannYhalfstep[u0, d1*hhy ] + fc[u0,v0] );
		v1= ADIXv@ (dirichletYhalfstep[v0, d2*hhy ] + gc[u0,v0] );
		{          ADIYu@ ( neumannXhalfstep[u1, d1*hhx ] + fc[u1,v1]),
		 setedges[ ADIYv@ (dirichletXhalfstep[v1, d2*hhx] + gc[u1,v1]), Last[bdata], a010] }];

RDstepDN[{u0_, v0_}, {fc_,gc_}, {d1_Real,d2_Real}, {hhx_,hhy_}]:= 
	Module[{u1,v1},
		u1= ADIXu@ (dirichletYhalfstep[u0, d1*hhy] + fc[u0,v0] );
		v1= ADIXv@ (  neumannYhalfstep[v0, d2*hhy] + gc[u0,v0] );
		{setedges[ ADIYu@ (dirichletXhalfstep[u1, d1*hhx ] + fc[u1,v1]), First[bdata], a010],
		           ADIYv@ (  neumannXhalfstep[v1, d2*hhx ] + gc[u1,v1])} ];
		 
RDstepNNN[{u0_, v0_, w0_}, {fc_,gc_,hc_}, {d1_,d2_,d3_}, {hhx_,hhy_}]:= 
	Module[{u1,v1,w1},
		u1= ADIXu@ (neumannYhalfstep[u0, d1*hhy] + fc[u0,v0,w0] );
		v1= ADIXv@ (neumannYhalfstep[v0, d2*hhy] + gc[u0,v0,w0] );
		w1= ADIXw@ (neumannYhalfstep[w0, d3*hhy] + hc[u0,v0,w0] );
		{   ADIYu@ (neumannXhalfstep[u1, d1*hhx] + fc[u1,v1,w1] ),
		    ADIYv@ (neumannXhalfstep[v1, d2*hhx] + gc[u1,v1,w1] ),
		    ADIYw@ (neumannXhalfstep[w1, d3*hhx] + hc[u1,v1,w1] ) }
		 ];
		 
RDstepNND[{u0_, v0_, w0_}, {fc_,gc_,hc_}, {d1_,d2_,d3_}, {hhx_,hhy_}]:= 
	Module[{u1,v1,w1},
		u1= ADIXu@ (  neumannYhalfstep[u0, d1*hhy ] + fc[u0,v0,w0] );
		v1= ADIXv@ (  neumannYhalfstep[v0, d2*hhy ] + gc[u0,v0,w0] );
		w1= ADIXw@ (dirichletYhalfstep[w0, d3*hhy ] + hc[u0,v0,w0] );
		{			ADIYu@ (  neumannXhalfstep[u1, d1*hhx ] + fc[u1,v1,w1] ),
		 			ADIYv@ (  neumannXhalfstep[v1, d2*hhx ] + gc[u1,v1,w1] ),
		  setedges[ ADIYw@ (dirichletXhalfstep[w1, d3*hhx ] + hc[u1,v1,w1] ), Last[bdata], a010]}
		 ];

RDstepNDD[{u0_, v0_, w0_}, {fc_,gc_,hc_}, {d1_,d2_,d3_}, {hhx_,hhy_}]:= 
	Module[{u1,v1,w1},
		u1= ADIXu@ (  neumannYhalfstep[u0, d1*hhy ] + fc[u0,v0,w0] );
		v1= ADIXv@ (dirichletYhalfstep[v0, d2*hhy ] + gc[u0,v0,w0] );
		w1= ADIXw@ (dirichletYhalfstep[w0, d3*hhy ] + hc[u0,v0,w0] );
		{          ADIYu@ (neumannXhalfstep[  u1, d1*hhx ] + fc[u1,v1,w1] ),
		 setedges[ ADIYv@ (dirichletXhalfstep[v1, d2*hhx ] + gc[u1,v1,w1] ), bdata[[2]],  a010],
		 setedges[ ADIYw@ (dirichletXhalfstep[w1, d3*hhx ] + hc[u1,v1,w1] ), Last[bdata], a010]}
		 ];
		 
RDstepDND[{u0_, v0_, w0_}, {fc_,gc_,hc_}, {d1_,d2_,d3_}, {hhx_,hhy_}]:= 
	Module[{u1,v1,w1},
		u1= ADIXu@ (dirichletYhalfstep[u0, d1*hhy ] + fc[u0,v0,w0]);
		v1= ADIXv@ (  neumannYhalfstep[v0, d2*hhy ] + gc[u0,v0,w0]);
		w1= ADIXw@ (dirichletYhalfstep[w0, d3*hhy ] + hc[u0,v0,w0]);
		{setedges[ ADIYu@ (dirichletXhalfstep[u1, d1*hhx ] + fc[u1,v1,w1]), First[bdata], a010],
		           ADIYv@ (  neumannXhalfstep[v1, d2*hhx ] + gc[u1,v1,w1]),
		 setedges[ ADIYw@ (dirichletXhalfstep[w1, d3*hhx ] + hc[u1,v1,w1]), Last[bdata],  a010]}
		 ];

RDstepDDD[{u0_, v0_, w0_}, {fc_,gc_,hc_}, {d1_,d2_,d3_}, {hhx_,hhy_}]:= 
	Module[{u1,v1,w1},
		u1= ADIXu@ (dirichletYhalfstep[u0, d1*hhy ] + fc[u0,v0,w0]);
		v1= ADIXv@ (dirichletYhalfstep[v0, d2*hhy ] + gc[u0,v0,w0]);
		w1= ADIXw@ (dirichletYhalfstep[w0, d3*hhy ] + hc[u0,v0,w0]);
		{setedges[ ADIYu@ (dirichletXhalfstep[u1, d1*hhx ] + fc[u1,v1,w1]), First[bdata], a010],
		 setedges[ ADIYv@ (dirichletXhalfstep[v1, d2*hhx ] + gc[u1,v1,w1]), bdata[[2]],   a010],
		 setedges[ ADIYw@ (dirichletXhalfstep[w1, d3*hhx ] + hc[u1,v1,w1]), Last[bdata],  a010] }
		 ];

RDstepDDN[{u0_, v0_, w0_}, {fc_,gc_,hc_}, {d1_,d2_,d3_}, {hhx_,hhy_}]:= 
	Module[{u1,v1,w1},
		u1= ADIXu@ (dirichletYhalfstep[u0, d1*hhy ] + fc[u0,v0,w0]);
		v1= ADIXv@ (dirichletYhalfstep[v0, d2*hhy ] + gc[u0,v0,w0]);
		w1= ADIXw@ (  neumannYhalfstep[w0, d3*hhy ] + hc[u0,v0,w0]);
		{setedges[ ADIYu@ (dirichletXhalfstep[u1, d1*hhx ] + fc[u1,v1,w1]), First[bdata], a010],
		 setedges[ ADIYv@ (dirichletXhalfstep[v1, d2*hhx ] + gc[u1,v1,w1]), bdata[[2]],   a010],
		           ADIYw@ (  neumannXhalfstep[w1, d3*hhx ] + hc[u1,v1,w1])}
		 ];

RDstepDNN[{u0_, v0_, w0_}, {fc_,gc_,hc_}, {d1_,d2_,d3_}, {hhx_,hhy_}]:= 
	Module[{u1,v1,w1},
		u1= ADIXu@ (dirichletYhalfstep[u0, d1*hhy ] + fc[u0,v0,w0]);
		v1= ADIXv@ (  neumannYhalfstep[v0, d2*hhy ] + gc[u0,v0,w0]);
		w1= ADIXw@ (  neumannYhalfstep[w0, d3*hhy ] + hc[u0,v0,w0]);
		{setedges[ ADIYu@ (dirichletXhalfstep[u1, d1*hhx ] + fc[u1,v1,w1]), First[bdata], a010],
		 		  ADIYv@ (  neumannXhalfstep[v1, d2*hhx ] + gc[u1,v1,w1]),
				   ADIYw@ (  neumannXhalfstep[w1, d3*hhx ] + hc[u1,v1,w1])}
		 ];
		 
RDstepNDN[{u0_, v0_, w0_}, {fc_,gc_,hc_}, {d1_,d2_,d3_}, {hhx_,hhy_}]:= 
	Module[{u1,v1,w1},
		u1= ADIXu@ (  neumannYhalfstep[u0, d1*hhy ] + fc[u0,v0,w0]);
		v1= ADIXv@ (dirichletYhalfstep[v0, d2*hhy ] + gc[u0,v0,w0]);
		w1= ADIXw@ (  neumannYhalfstep[w0, d3*hhy ] + hc[u0,v0,w0]);
		{		 ADIYu@ (  neumannXhalfstep[u1, d1*hhx ] + fc[u1,v1,w1]),
		 setedges[ADIYv@ (dirichletXhalfstep[v1, d2*hhx ] + gc[u1,v1,w1]), bdata[[2]], a010],
				  ADIYw@ (  neumannXhalfstep[w1, d3*hhx ] + hc[u1,v1,w1])}
		 ];
		 
RDstepPeriodic3[{u0_, v0_, w0_}, {fc_,gc_,hc_}, {d1_,d2_,d3_}, {hhx_,hhy_}]:= 
	Module[{u1,v1,w1},
		u1= ADIXu@(periodicYhalfstep[u0, d1*hhy ] + fc[u0,v0,w0] );
		v1= ADIXv@(periodicYhalfstep[v0, d2*hhy ] + gc[u0,v0,w0] );
		w1= ADIXw@(periodicYhalfstep[w0, d3*hhy ] + hc[u0,v0,w0] );
	 	{  ADIYu@(periodicXhalfstep[u1, d1*hhx ] + fc[u1,v1,w1] ),
		    ADIYv@(periodicXhalfstep[v1, d2*hhx ] + gc[u1,v1,w1] ),
		    ADIYw@(periodicXhalfstep[w1, d3*hhx ] + hc[u1,v1,w1] ) }
		 ];
		 		
BlockDiag[(mat_)?MatrixQ.., k_] :=
  SparseArray[Flatten[MapThread[Array[List, #1, #2 + 1]&, ({#1, Most[FoldList[Plus, {0, 0}, #1]]}&)[
       PadRight[{{k,k}}, m*n/k, {{k,k}} ] ] ], 2] -> PadRight[Flatten[mat], m*n*k, Flatten[mat]]];

sparray[bc_,d_,h_,k_]:=
    Which[
        bc==="Periodic", BlockDiag[
        	Compile[{{kk, _Integer}},
        		Table[ Which[i==j, 1+d*h, Abs[i-j]==1, -.5d*h, Abs[i-j]==k-1, -.5d*h, True, 0],
        	  	{i, kk}, {j, kk}] ][k], k
        	  	],
		bc==="Dirichlet", BlockDiag[
        	Compile[{{kk, _Integer}},
        		Table[ Which[i==j==1, 1, i==j==k, 1, i==j, 1+d*h, Abs[i-j]==1, -.5d*h, True, 0],
        	  	{i, kk}, {j, kk}] ][k], k
        	  	],
		bc==="Neumann", BlockDiag[
        	Compile[{{kk, _Integer}},
        		Table[ Which[i==j==1, 1+.5d*h, i==j==k, 1+.5d*h, i==j, 1+d*h, Abs[i-j]==1, -.5d*h, True, 0],
        	  	{i, kk}, {j, kk}] ][k], k
        	  	]
        ];

makeLSolvers[{bcu_,bcv_},{du_,dv_},h_]:=
        LinearSolve/@
            { \[Tau]++; sparray[bcu, du, First@h, n],
			  \[Tau]++; sparray[bcu, du, Last@h, m],
			  \[Tau]++; sparray[bcv, dv, First@h, n],
			  \[Tau]++; sparray[bcv, dv, Last@h, m] };
            
makeLSolvers[{bcu_,bcv_,bcw_},{du_,dv_,dw_},h_]:=
        LinearSolve/@
            { \[Tau]++; sparray[bcu, du, First@h, n],
			  \[Tau]++; sparray[bcu, du, Last@h, m],
			  \[Tau]++; sparray[bcv, dv, First@h, n],
			  \[Tau]++; sparray[bcv, dv, Last@h, m],
			  \[Tau]++; sparray[bcw, dw, First@h, n],
			  \[Tau]++; sparray[bcw, dw, Last@h, m] };

fasterDensityPlot[u_?MatrixQ,opts___]:=
	Module[{pr, m=Min[u], M=Max[u]+$MachineEpsilon}, 
		pr= (PlotRange/.Flatten[{opts}])/.All->{m-.01Abs[m], M+.01Abs[M]};
			Graphics[
			Raster[u,{{0,0},{1,1}}, pr,
				Sequence[FilterRules[{opts},Options[Raster]]] ], 
			Sequence[DeleteCases[FilterRules[{opts},Options[Graphics]],(PlotRange->_)|(PlotRangePadding->_)]]
			 ] ];

RDDensityPlots[{u_?MatrixQ,v_?MatrixQ}/;Dimensions[u]===Dimensions[v], opts___Rule ]:= 
		Module[{x,y}, RDDensityPlots[{1,2},{x,y},{u,v},{1,1},{0,1,1}, NPlots->0, opts]];
RDDensityPlots[{u_?MatrixQ,v_?MatrixQ,w_?MatrixQ}/;Dimensions[u]===Dimensions[v]===Dimensions[w], opts___Rule ]:= 
		Module[{x,y,z}, RDDensityPlots[{1,2,3},{x,y,z},{u,v,w},{1,1,1},{0,1,1}, NPlots->0, opts]]

		
(******* Definition of RDDensityPlots for 2-component systems *********************************)

RDDensityPlots[ {f_,g_}, {u_Symbol,v_Symbol}, 
			{u0_?MatrixQ, v0_?MatrixQ}/; Dimensions[u0]===Dimensions[v0],
			{d1_?NumericQ, d2_?NumericQ}/; And@@NonNegative/@{d1,d2},
			{t0_?NumericQ, t1_?NumericQ, dt_?NumericQ}/; And[t0 \[Element] Reals, t1 \[Element] Reals, Positive[dt]], 
			opts___Rule] :=
	
	Module[{a,a101,ar,b,bc,bdf,bnds,c,d,fc,gc,grsp,hh,ims,kk,lbls,ldpopts,maxu,maxv,mf1,mf2,
				minu,minv,mMuv,nplots,pfg,plotter,pr,step,uvals,vvals,wait,lb,ub,\[Delta]t,ws},

	 (* Not local to this module: m, n, \[Tau], a010, bdata, uLsolvex, uLsolvey, vLsolvex, vLsolvey *)

		If[Not[NumericQ[ f^2+g^2 /. {u->RandomReal[],v->RandomReal[]}]],
				Message[RDDensityPlots::nnfns]; Return[$Failed]];
		If[Not[And@@(MatchQ[N[#],_Real]&)/@Flatten[{u0,v0}]],
				Message[RDDensityPlots::nndata]; Return[$Failed]];

		{m,n} = Dimensions[u0]; 

		pfg= PerformanceGoal/.{opts}/.Options[RDDensityPlots];
		plotter= If[pfg==="Speed", fasterDensityPlot, ListDensityPlot];

		nplots= NPlots/.{opts}/.Options[RDDensityPlots];
		If[nplots==Automatic, nplots=If[pfg==="Speed", Clip[Abs[t1-t0], {10,100}], Clip[Abs[t1-t0], {5,25}] ] ];
		If[NumericQ[nplots]&& nplots>=.5, nplots=Round[nplots]];
		If[nplots!=0, nplots=Max[nplots,1]];
		If[nplots>0, 
				wait= Max[Round[Abs[t1-t0]/nplots/dt],1]; \[Delta]t= N[Abs[t1-t0]/nplots/wait]
			];  
  \[Tau]=0; 
	mf1= If[m*n > 3000 && nplots>0, Monitor, #1&];
	mf1[
		ldpopts= FilterRules[Join[{opts}, DeleteCases[Options[RDDensityPlots], AspectRatio->_ ] ],
							 Options[ListDensityPlot] ];
		grsp= Spacings/.{opts}/. Options[RDDensityPlots]/. Automatic->0;
		If[Not[NumericQ[grsp]]||Not[Positive[grsp]], grsp=0 ];  
		lbls = PlotLabel/.{opts}/.Options[RDDensityPlots ];
		If[ Not[MatchQ[lbls,{_,_}]], lbls={None,None} ];
		
		bc = Flatten[{ BoundaryConditions/.{opts}/.Options[RDDensityPlots] } ]; 
		bc= If[Length[bc]===1, Flatten[{bc,bc}], bc];
		bc = Replace[bc, x_String:>(x->(0&)),1];
		bdf = bc/.(_->val_):>val;
		bc = bc/.(o_->_):>o; 

		ws= First[AbsoluteCurrentValue[WindowSize]]; 
		ims= ImageSize/.{opts}/.Options[RDDensityPlots]/.{Tiny:>.25ws,Small:>.33ws,Medium:>.5ws,Large:>.75ws};
		ims= If[NumericQ[ims]&&ims>40,ims,400];

		bnds= Bounds/.{opts}/.Options[RDDensityPlots];
		{lb,ub}= Which[NumericQ[bnds]&&Positive[bnds], {0,bnds},
						MatchQ[bnds, {x_?NumericQ,y_?NumericQ}/;And@@(# \[Element] Reals&/@{x,y})], bnds,
						True, {0,1000}  ];

		minu = Min[u0]; maxu = Max[u0];  
		minv = Min[v0]; maxv = Max[v0];
		mMuv={ {minu-.01Abs[minu], maxu+.01Abs[maxu]},
				{minv-.01Abs[minv], maxv+.01Abs[maxv]}  };

		pr = PlotRange/.{opts}/.Options[RDDensityPlots];
		pr = If[ pr===Automatic||pr===All||MatchQ[pr,{_?NumericQ,_?NumericQ}], {pr, pr}, pr  ];  
		If[ Not[MatchQ[pr, {({_?NumericQ,_?NumericQ}|Automatic|All) ..}]], 
			Message[RDDensityPlots::pr]; Abort[]];
		pr= MapIndexed[If[#1===Automatic, mMuv[[Sequence@@#2]], #1]&, pr];
		
		{{a,b},{c,d}}= N[ Domain/.{opts}/. Options[RDDensityPlots]/.
							{{x1_,Automatic},{y1_,y2_}}:>{{x1,(n-1)/(m-1)},{y1,y2}} ];

		If[Not[MatchQ[{{a,b},{c,d}},{{_?NumericQ ..} ..}]], {{a,b},{c,d}} = {{0,(n-1)/(m-1)},{0,1}}];
		
		ar = (d-c)/(b-a);
		
		hh = \[Delta]t/{((b-a)/(n-1))^2,((d-c)/(m-1))^2};
	
		{fc,gc}= Compile[{{u,_Real,2},{v,_Real,2}}, .5*\[Delta]t*#]&/@{f,g}; 
		
		a010 = Table[ If[ i==1||j==1||i==m||j==n, 0, 1 ], {i,m}, {j,n}];
		a101 = Mod[1+a010,2];
		
		step= Which[
				bc==={"Periodic","Periodic"},   RDstepPeriodic,
				bc==={"Dirichlet","Dirichlet"}, RDstepDD,
				bc==={"Neumann","Neumann"},     RDstepNN,
				bc==={"Dirichlet","Neumann"},   RDstepDN,
				bc==={"Neumann","Dirichlet"},   RDstepND,
				True, Message[RDDensityPlots::bc]; Return[$Failed] ];
						
		If[ MemberQ[bc,"Dirichlet"],
				bdata= { If[bc[[1]]==="Dirichlet", 
							a101*Table[First[bdf][x,y], {y,c,d,(d-c)/(m-1)}, {x,a,b,(b-a)/(n-1)} ], {} ], 
						 If[bc[[2]]==="Dirichlet", 
							a101*Table[ Last[bdf][x,y], {y,c,d,(d-c)/(m-1)}, {x,a,b,(b-a)/(n-1)} ], {} ] }
			]; 

		If[nplots>0, {uLsolvex,uLsolvey,vLsolvex,vLsolvey} = makeLSolvers[bc, {d1,d2}, hh]];
  \[Tau]++;
		{uvals,vvals}=
			N[Which[
			bc==={"Dirichlet","Dirichlet"}, {setedges[u0,First[bdata],a010], setedges[v0,Last[bdata],a010]},
			bc==={"Neumann","Neumann"},     {u0,v0},
			bc==={"Periodic","Periodic"},   {periodicize[u0], periodicize[v0]},
			bc==={"Dirichlet","Neumann"},   {setedges[u0, First[bdata], a010], v0},
			bc==={"Neumann","Dirichlet"},   {u0, setedges[v0, Last[bdata], a010]}
			] ], 

		 Row[{Spacer[20], Style["    setting up linear solvers", FontFamily->"Trebuchet MS", 10], "   ",
			ProgressIndicator[\[Tau]/5]}] 
				];	

(**)  DynamicModule[{\[Omega]=wait, t=t0-\[Delta]t, j=0, s=ims, \[Delta]\[Tau]=\[Delta]t, arat=ar, pic, minmax, RowOrColumn, row, col, asm, sp=grsp},	
			mf2=If[(nplots*wait<40 && pfg==="Speed")||nplots==0, #1&, Monitor];
			row=(ar>.6);col=(ar<=.6);
			RowOrColumn= If[col, Column, Row]; 
			kk=0; 
			t= t+\[Delta]\[Tau];
			{uvals,vvals}={u0,v0};
			minmax[0]={Min[#],Max[#]}&/@{uvals,vvals};
			pic[0]= Table[ plotter[{uvals,vvals}[[i]], PlotRange->pr[[i]], 
								PlotLabel->lbls[[i]], Sequence@@ldpopts, AspectRatio->ar ],
						{i,1,2}];
			If[nplots>0,	
				mf2[ Do[
						{uvals,vvals}= 
							Nest[
								Module[{test=(lb<=Min[#]<=Max[#]<=ub&), next=step[t=t+\[Delta]\[Tau];kk++; #, {fc,gc}, {d1,d2}, hh]},
									If[And@@(test/@next), 
										Chop[next,$MachineEpsilon], 
										Message[RDDensityPlots::outbnds]; Abort[] ] ]&,
									{uvals,vvals}, \[Omega]];
							minmax[k]={Min[#],Max[#]}&/@{uvals,vvals};
							pic[k]= Table[plotter[{uvals,vvals}[[i]], 
										PlotRange->pr[[i]], PlotLabel->lbls[[i]], Sequence@@ldpopts, AspectRatio->ar ],
									{i,1,2}],
					{k,1,nplots} ],
				
					Row[{Spacer[20], Style["computing solution", FontFamily->"Trebuchet MS",10], "   ",
					ProgressIndicator[kk/(nplots wait+.1)], "   ", 
					Style["\[ScriptT] = ", FontFamily->"Trebuchet MS", 11],
					Style[TraditionalForm[t0+kk*\[Delta]\[Tau]], FontFamily->"Trebuchet MS",11]}] 
					
				] ]; 

			RDL`LastValues={uvals,vvals}; 
			RDL`Images=Table[pic[i],{i,0,nplots}];
			asm=DefaultAnimationSpeedMultiplier/.{opts}/.Options[RDDensityPlots];
			If[Not[NumericQ[asm]],asm=1.5]; asm= Clip[asm,{.01,100}];
			Column[{
				If[nplots>0,
						Row[{Spacer[10], Animator[Dynamic[j], {0,nplots,1}, asm/Min[1,(\[Delta]\[Tau]*\[Omega])], 
									Appearance->Small, AnimationRunning->False, AnimationRepetitions->1 ],
							"   ", Style["\[ScriptT] = ", FontFamily->"Trebuchet MS", 11],
							Style[TraditionalForm[Dynamic[ t0+j*\[Omega]*\[Delta]\[Tau] ]], FontFamily->"Trebuchet MS", 11] } ], 
						Sequence@@{}], 
				Dynamic@RowOrColumn[
						Table[Column[{Show[pic[j][[i]], ImageSize->If[col, s, s/2.01]],
								Row[{Style[PaddedForm[Round[ minmax[j][[i]], .00001],{6,4}], 10, FontFamily->"Arial"]},
										 ImageSize-> s/2.01, Alignment->Center]
									}, Spacings->{If[col,0,sp],0}, Frame->True, FrameStyle->Opacity[0], Alignment->Center],
							{i,1,2}],
								If[col, Spacings->{0,sp}, Sequence@@{}] 
							] }, Spacings-> If[col,sp,.2+sp/5]]  
			]	
		]

(************* Definition of RDDensityPlots for 3-component systems *******************************)
		
RDDensityPlots[ {f_,g_,h_}, {u_Symbol,v_Symbol,w_Symbol},
		{u0_?MatrixQ, v0_?MatrixQ, w0_?MatrixQ}/; Dimensions[u0]===Dimensions[v0]===Dimensions[w0],
		{d1_?NumericQ, d2_?NumericQ, d3_?NumericQ}/; And@@NonNegative/@{d1,d2,d3}, 
		{t0_?NumericQ, t1_?NumericQ, dt_?NumericQ}/; And[t0 \[Element] Reals, t1 \[Element] Reals, Positive[dt]], 
		opts___Rule ]:=
		
	Module[{a,a101,ar,b,bc,bdf,bnds,c,d,fc,gc,grsp,hc,hh,ims,kk,lbls,ldpopts,maxu,maxv,maxw,mf1,mf2,
			minu,minv,minw,mMuvw,nplots,pfg,plotter,pr,step,uvals,vvals,wait,wvals,lb,ub,\[Delta]t,ws},

		(* Not local to this module: m, n, \[Tau], a010, bdata, uLsolvex, uLsolvey, vLsolvex, vLsolvey *)

		If[Not[NumericQ[f^2+g^2+h^2/.{u->RandomReal[],v->RandomReal[],w->RandomReal[]}]], 
				Message[RDDensityPlots::nnfns]; Return[$Failed]];
		If[Not[And@@(MatchQ[N[#],_Real]&)/@Flatten[{u0,v0,w0}]],
				Message[RDDensityPlots::nndata]; Return[$Failed]];
   \[Tau]=0;
		{m,n} = Dimensions[u0];

		pfg= PerformanceGoal/.{opts}/.Options[RDDensityPlots];
		plotter= If[pfg==="Speed", fasterDensityPlot, ListDensityPlot];

		nplots= NPlots/.{opts}/.Options[RDDensityPlots];
		If[nplots==Automatic, nplots= If[pfg==="Speed", Clip[Abs[t1-t0], {10,100}], Clip[Abs[t1-t0], {5,25}] ] ];
		If[NumericQ[nplots] && nplots>=.5, nplots= Round[nplots] ];
		If[nplots!=0, nplots=Max[nplots,1]];
		If[nplots>0, 
				wait= Max[Round[Abs[t1-t0]/nplots/dt],1]; \[Delta]t= N[Abs[t1-t0]/nplots/wait]
			]; 

		mf1= If[m*n > 2500 && nplots>0, Monitor, #1&];
	mf1[ 		
		ldpopts= FilterRules[Join[{opts}, 
							DeleteCases[Options[RDDensityPlots], AspectRatio->_ ]], 
							Options[ListDensityPlot]];
		grsp= Spacings/.{opts}/. Options[RDDensityPlots]/. Automatic->0;
		If[Not[NumericQ[grsp]]||Not[Positive[grsp]], grsp=0 ]; 
		grsp= Spacings/.{opts}/.Options[RDDensityPlots];
		lbls = PlotLabel/.{opts}/.Options[RDDensityPlots];
		If[ Not[MatchQ[lbls,{_,_,_}]], lbls={None,None,None}];

		bnds= Bounds/.{opts}/.Options[RDDensityPlots];
		{lb,ub}= Which[NumericQ[bnds]&&Positive[bnds], {0,bnds},
						MatchQ[bnds, {x_?NumericQ,y_?NumericQ}/;And@@(# \[Element] Reals&/@{x,y})], Sort[bnds],
						True, {0,1000}  ];
				
		bc = Flatten[{BoundaryConditions/.{opts}/.Options[RDDensityPlots]}];  
		bc = If[Length[bc]===1, Flatten@{bc,bc,bc}, bc];
		bc = Replace[bc, x_String:>(x->(0&)),1];
		bdf= bc/.(_->val_):>val;
		bc = bc/.(o_->_):>o;

		ws = First[CurrentValue[WindowSize]]; 
		ims= ImageSize/.{opts}/.Options[RDDensityPlots]/.{Tiny:>.25ws,Small:>.33ws,Medium:>.5ws,Large:>.75ws};
		ims= If[NumericQ[ims]&&ims>40,ims,400];

		minu=Min[u0]; maxu=Max[u0];
		minv=Min[v0]; maxv=Max[v0];
		minw=Min[w0]; maxw=Max[w0]; 
		mMuvw={ {minu-.01Abs[minu], maxu+.01Abs[maxu]},
				{minv-.01Abs[minv], maxv+.01Abs[maxv]},
				{minw-.01Abs[minw], maxw+.01Abs[maxw]}  };

		pr= PlotRange/.{opts}/.Options[RDDensityPlots];
		pr= If[pr===Automatic||pr===All||VectorQ[pr], {pr,pr,pr}, pr]; 
		If[Not[MatchQ[pr, {({_?NumericQ,_?NumericQ}|Automatic|All) ..}]],
			 Message[RDDensityPlots::pr]; Abort[]];
		pr= MapIndexed[If[#1===Automatic, mMuvw[[Sequence@@#2]], #1]&, pr];


		{{a,b},{c,d}} = N[ Domain/.{opts}/.Options[RDDensityPlots]/.
									{{x1_,Automatic},{y1_,y2_}} :> {{x1,(n-1)/(m-1)},{y1,y2}} ];
		
		If[Not[MatchQ[{{a,b},{c,d}},{{_?NumericQ ..} ..}]], {{a,b},{c,d}} = {{0,(n-1)/(m-1)},{0,1}}];
		
		ar = (d-c)/(b-a);
		
		hh = \[Delta]t/{((b-a)/(n-1))^2,((d-c)/(m-1))^2}; 
  \[Tau]++ ;	
		{fc,gc,hc}= Compile[{{u,_Real,2},{v,_Real,2},{w,_Real,2}}, .5*\[Delta]t*#]&/@{f,g,h}; 	
		
		a010 = Table[ If[ i==1||j==1||i==m||j==n, 0, 1], {i,m}, {j,n}]; 
		a101 = Mod[1+a010,2]; 
  \[Tau]++ ;	
		step= Which[
			bc==={"Dirichlet","Dirichlet","Dirichlet"}, RDstepDDD,
			bc==={"Neumann","Dirichlet","Dirichlet"},   RDstepNDD,
			bc==={"Dirichlet","Neumann","Dirichlet"},   RDstepDND,
			bc==={"Neumann","Neumann","Dirichlet"},     RDstepNND,
			bc==={"Dirichlet","Dirichlet","Neumann"},   RDstepDDN,
			bc==={"Neumann","Dirichlet","Neumann"},     RDstepNDN,
			bc==={"Dirichlet","Neumann","Neumann"},     RDstepDNN,
			bc==={"Neumann","Neumann","Neumann"},       RDstepNNN,
			bc==={"Periodic","Periodic","Periodic"},    RDstepPeriodic3,
			True, Message[RDDensityPlots::bc]; Return[$Failed]
			];	
		
		If[MemberQ[bc, "Dirichlet"],
			bdata= {If[bc[[1]]==="Dirichlet",
						a101*Table[First[bdf][x,y], {y,c,d,(d-c)/(m-1)}, {x,a,b,(b-a)/(n-1)} ], {} ],
					If[bc[[2]]==="Dirichlet", 
						a101*Table[  bdf[[2]][x,y], {y,c,d,(d-c)/(m-1)}, {x,a,b,(b-a)/(n-1)} ], {} ],
					If[bc[[3]]==="Dirichlet", 
						a101*Table[ Last[bdf][x,y], {y,c,d,(d-c)/(m-1)}, {x,a,b,(b-a)/(n-1)} ], {} ] 
				} ];

		If[nplots>0, 
			{uLsolvex,uLsolvey,vLsolvex,vLsolvey,wLsolvex,wLsolvey} = makeLSolvers[bc, {d1,d2,d3}, hh] ]; 
					
		{uvals,vvals,wvals}= N[
			Which[
			bc==={"Dirichlet","Dirichlet","Dirichlet"},
				{setedges[u0,bdata[[1]],a010], setedges[v0,bdata[[2]],a010], setedges[w0,bdata[[3]],a010]},
			bc==={"Neumann","Dirichlet","Dirichlet"},
				{u0, setedges[v0,bdata[[2]],a010], setedges[w0,bdata[[3]],a010]},
			bc==={"Dirichlet","Neumann","Dirichlet"},
				{setedges[u0,bdata[[1]],a010], v0, setedges[w0,bdata[[3]],a010]},
			bc==={"Neumann","Neumann","Dirichlet"},
				{u0, v0, setedges[w0,bdata[[3]],a010]},
			bc==={"Neumann","Dirichlet","Neumann"},
				{u0, setedges[v0,bdata[[2]], a010], w0},
			bc==={"Dirichlet","Dirichlet","Neumann"},
				{setedges[u0,bdata[[1]],a010], setedges[v0,bdata[[2]],a010], w0},
			bc==={"Dirichlet","Neumann","Neumann"},
				{setedges[u0,bdata[[1]],a010], v0, w0},
			bc==={"Neumann","Neumann","Neumann"},
				{u0, v0, w0},
			bc==={"Periodic","Periodic","Periodic"},
				{periodicize[u0], periodicize[v0], periodicize[w0]}
			] ] ; \[Tau]++ , 

		 Row[{Spacer[20], Style["    setting up stuff", FontFamily->"Trebuchet MS",10], "   ",
			ProgressIndicator[\[Tau]/9]}] 
		];

(***)  DynamicModule[{\[Omega]=wait, j=0, \[Delta]\[Tau]=\[Delta]t, s=ims, arat=ar, pic, minmax, RowOrColumn, row,col, asm, sp=grsp},	
			mf2=If[(nplots*wait<40 && pfg==="Speed")||nplots==0, #1&, Monitor];
			row=(ar>.4);col=(ar<=.4);
			RowOrColumn= If[col, Column, Row]; 
			kk=0; 
			{uvals,vvals,wvals}= {u0,v0,w0};
			minmax[0]= {Min[#],Max[#]}&/@{uvals,vvals,wvals};
			pic[0]= Table[ plotter[{uvals,vvals,wvals}[[i]], PlotRange->pr[[i]], 
								PlotLabel->lbls[[i]], Sequence@@ldpopts, AspectRatio->ar ],
						{i,1,3}];
			If[nplots>0,
				mf2[ Do[
						{uvals,vvals,wvals}= 
							Nest[
								Module[{test=(lb<=Min[#]<=Max[#]<=ub&), next=step[kk++; #, {fc,gc,hc}, {d1,d2,d3}, hh]},
									If[And@@(test/@next), 
										Chop[next,$MachineEpsilon],
										Message[RDDensityPlots::outbnds]; Abort[] ] ]&,
									{uvals,vvals,wvals}, \[Omega]];
							minmax[k]= {Min[#],Max[#]}&/@{uvals,vvals,wvals};
							pic[k]= Table[plotter[{uvals,vvals,wvals}[[i]], 
										PlotRange->pr[[i]], PlotLabel->lbls[[i]], Sequence@@ldpopts, AspectRatio->ar ],
									{i,1,3}],
					{k,1,nplots} ],				
				
					Row[{Spacer[20], Style["computing solution", FontFamily->"Trebuchet MS", 10 ], "   ",
					ProgressIndicator[kk/(nplots wait+.1)], "   ", 
					Style["\[ScriptT] = ", FontFamily->"Trebuchet MS", 11],
					Style[TraditionalForm[t0+kk*\[Delta]t], FontFamily->"Trebuchet MS", 11]}] 
					
				] ];

			RDL`LastValues= {uvals,vvals,wvals}; 
			RDL`Images= Table[pic[i],{i,0,nplots}];
			asm=DefaultAnimationSpeedMultiplier/.{opts}/.Options[RDDensityPlots];
			If[Not[NumericQ[asm]],asm=1.5]; asm= Clip[asm,{.01,100}];

			Column[{
				If[nplots>0,
						Row[{Spacer[10], Animator[Dynamic[j], {0,nplots,1}, asm/Min[1,\[Delta]t*\[Omega]], 
										Appearance->Small, AnimationRunning->False, AnimationRepetitions->1],
							 "   ", Style["\[ScriptT] = ", FontFamily->"Trebuchet MS", 11],
							Style[TraditionalForm[Dynamic[ t0+j*\[Omega]*\[Delta]\[Tau] ]], FontFamily->"Trebuchet MS", 11 ] } ], 
						Sequence@@{}], 
				Dynamic@RowOrColumn[
						Table[Column[{Show[pic[j][[i]], ImageSize->If[col, s, s/3.01] ],
								Row[{Style[PaddedForm[Round[ minmax[j][[i]], .00001], {6, 4} ], 10, FontFamily->"Arial"] },
										 ImageSize-> s/3.01, Alignment->Center]
									}, Spacings->{If[col,0,sp],0}, Frame->True, FrameStyle->Opacity[0], Alignment->Center],
						{i, 1, 3}]
							] }, Spacings->If[col,sp,.2+sp/5]]  
			]	
		]

End[];

Protect[RDDensityPlots];

EndPackage[];

Print["                    ok, got it!"]

