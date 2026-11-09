Group{
  DefineGroup[
    DomainS0_Mag // Stranded model if imposing source with js0[]
  ];
}

FunctionSpace {

  { Name Hcurl_a_Mag_2D; Type Form1P;
    BasisFunction {
      { Name se; NameOfCoef ae; Function BF_PerpendicularEdge;
        Support Domain_Mag; Entity NodesOf[ All ]; }
      If (Flag_Degree_a == 2)
        { Name se2; NameOfCoef ae2; Function BF_PerpendicularEdge_2E;
          Support Domain_Mag; Entity EdgesOf[ All ]; }
      EndIf
    }
    Constraint {
      { NameOfCoef ae;
        EntityType NodesOf; NameOfConstraint MagneticVectorPotential_2D; }
      If (Flag_Degree_a == 2)
	    { NameOfCoef ae2; // Only OK if homogeneous BC, otherwise specify zero-BC
          EntityType EdgesOf; NameOfConstraint MagneticVectorPotential_2D; }
      EndIf
    }
  }

  { Name Hregion_i_2D ; Type Vector ;
    BasisFunction {
      { Name sr ; NameOfCoef ir ; Function BF_RegionZ ;
        Support DomainS_Mag ; Entity DomainS_Mag ; }
    }
    GlobalQuantity {
      { Name Is ; Type AliasOf        ; NameOfCoef ir ; }
      { Name Us ; Type AssociatedWith ; NameOfCoef ir ; }
    }
    Constraint {
      { NameOfCoef Us ; EntityType Region ; NameOfConstraint Voltage_2D ; }
      { NameOfCoef Is ; EntityType Region ; NameOfConstraint Current_2D ; }
    }
  }

  // Gradient of Electric scalar potential (2D)
  { Name Hregion_u_Mag_2D ; Type Form1P ;
    BasisFunction {
      { Name sr ; NameOfCoef ur ; Function BF_RegionZ ;
        Support DomainC_Mag ; Entity DomainC_Mag ; }
    }
    GlobalQuantity {
      { Name U ; Type AliasOf        ; NameOfCoef ur ; }
      { Name I ; Type AssociatedWith ; NameOfCoef ur ; }
    }
    Constraint {
      { NameOfCoef U ;
        EntityType Region ; NameOfConstraint Voltage_2D ; }
      { NameOfCoef I ;
        EntityType Region ; NameOfConstraint Current_2D ; }
    }
  }


}

Formulation {

  { Name Magnetodynamic_a_2D; Type FemEquation; // Magnetodynamics
    Quantity {
      { Name a;  Type Local; NameOfSpace Hcurl_a_Mag_2D; }

      // massive conductors (source or not)
      { Name ur ; Type Local  ; NameOfSpace Hregion_u_Mag_2D ; }
      { Name I  ; Type Global ; NameOfSpace Hregion_u_Mag_2D [I] ; }
      { Name U  ; Type Global ; NameOfSpace Hregion_u_Mag_2D [U] ; }

       // stranded conductors (source)
      { Name ir ; Type Local  ; NameOfSpace Hregion_i_2D ; }
      { Name Us ; Type Global ; NameOfSpace Hregion_i_2D[Us] ; }
      { Name Is ; Type Global ; NameOfSpace Hregion_i_2D[Is] ; }
    }

    Equation {
      Galerkin { [ nu[] * Dof{d a} , {d a} ];
        In Domain_Mag; Jacobian Vol; Integration I1; }
      Galerkin { DtDof [ sigma[] * Dof{a} , {a} ];
        In DomainC_Mag; Jacobian Vol; Integration I1; }

      Galerkin { [ sigma[] * Dof{ur}, {a} ] ;
        In DomainC_Mag ; Jacobian Vol ; Integration I1 ; }
      Galerkin { DtDof [ sigma[] * Dof{a} , {ur} ] ;
        In DomainC_Mag ; Jacobian Vol ; Integration I1 ; }
      Galerkin { [ sigma[] * Dof{ur}, {ur} ] ;
        In DomainC_Mag ; Jacobian Vol ; Integration I1 ; }

      GlobalTerm { [ Dof{I} , {U} ] ; In DomainCWithI_Mag ; }

      // Either you impose directly the function js0[]
      Galerkin { [ -js0[] , {a} ];
        In DomainS0_Mag; Jacobian Vol; Integration I1; }

      // or you use the constraints => allows accounting for sigma[]
      Galerkin { [ -Ns[]/Sc[] * Dof{ir}, {a} ] ;
        In DomainS_Mag ; Jacobian Vol ; Integration I1 ; }
      Galerkin { DtDof [ Ns[]/Sc[] * Dof{a}, {ir} ] ;
        In DomainS_Mag ; Jacobian Vol ; Integration I1 ; }

      Galerkin { [ Ns[]/Sc[] / sigma[] * Ns[]/Sc[]* Dof{ir} , {ir} ] ; // resistance term
        In DomainS_Mag ; Jacobian Vol ; Integration I1 ; }
      //GlobalTerm { [ Rdc * Dof{Is} , {Is} ] ; In DomainS ; } // OR this resitance term
      GlobalTerm { [ Dof{Us}, {Is} ] ; In DomainS_Mag ; }

    }
  }

}


Resolution {

  { Name Magnetodynamics;
    System {
      { Name Sys_Mag; NameOfFormulation Magnetodynamic_a_2D;
        Type Complex; Frequency Freq; }
    }
    Operation {
      CreateDir["res"];

      InitSolution[Sys_Mag];
      Generate[Sys_Mag]; Solve[Sys_Mag]; SaveSolution[Sys_Mag];
      PostOperation[Mag_Maps];
      PostOperation[Mag_Global];
    }
  }

}


PostProcessing {

  { Name Magnetodynamic_a_2D; NameOfFormulation Magnetodynamic_a_2D;
    PostQuantity {
      { Name a; Value { Term { [ {a} ]; In Domain_Mag; Jacobian Vol; } } }
      { Name az; Value { Term { [ CompZ[{a}] ]; In Domain_Mag; Jacobian Vol; } } }
      { Name b; Value { Term { [ {d a} ]; In Domain_Mag; Jacobian Vol; } } }
      { Name bm; Value { Term { [ Norm[{d a}] ]; In Domain_Mag; Jacobian Vol; } } }

      { Name j; Value {
          Term { [ -sigma[]*(Dt[{a}]+{ur})]; In DomainC_Mag; Jacobian Vol; }
          Term { [ js0[] ]; In DomainS0_Mag; Jacobian Vol; }
          Term { [ Ns[]/Sc[]*{ir} ]; In DomainS_Mag; Jacobian Vol; }
        } }
      { Name jz; Value {
          Term { [ CompZ[-sigma[]*(Dt[{a}]+{ur})] ]; In DomainC_Mag; Jacobian Vol; }
          Term { [ CompZ[ js0[] ] ]; In DomainS0_Mag; Jacobian Vol; }
          Term { [ CompZ[ Ns[]/Sc[]*{ir} ]]; In DomainS_Mag; Jacobian Vol; }
        } }
      { Name jm; Value {
          Term { [ Norm[-sigma[]*(Dt[{a}]+{ur})] ]; In DomainC_Mag; Jacobian Vol; }
          Term { [ Norm[ js0[] ] ]; In DomainS0_Mag; Jacobian Vol; }
          Term { [ Norm[ Ns[]/Sc[]*{ir} ]]; In DomainS_Mag; Jacobian Vol; }
        } }

      { Name rhoj2; Value {
          Term { [ 0.5*sigma[]*SquNorm[Dt[{a}]+{ur}] ]; In DomainC_Mag; Jacobian Vol; }
          Term { [ 0.5/sigma[]*SquNorm[js0[]] ]; In DomainS0_Mag ; Jacobian Vol; }
          Term { [ 0.5/sigma[]*SquNorm[Ns[]/Sc[]*{ir}] ]; In DomainS_Mag; Jacobian Vol; }
        }
      }


      { Name JouleLosses; Value {
          Integral { [ 0.5*sigma[]*SquNorm[Dt[{a}]] ]   ; In DomainC_Mag  ; Jacobian Vol ; Integration I1 ; }
          Integral { [ 0.5/sigma[]*SquNorm[js0[]] ]          ; In DomainS0_Mag ; Jacobian Vol ; Integration I1 ; }
          Integral { [ 0.5/sigma[]*SquNorm[Ns[]/Sc[]*{ir}] ] ; In DomainS_Mag  ; Jacobian Vol ; Integration I1 ; }
        }
      }

      { Name U ; Value {
          Term { [ {U} ]  ; In DomainC_Mag ; }
          Term { [ {Us} ] ; In DomainS_Mag ; }
        }
      }

      { Name I ; Value {
          Term { [ {I} ]  ; In DomainC_Mag ; }
          Term { [ {Is} ] ; In DomainS_Mag ; }
        }
      }

      { Name S ; Value {
          Term { [ {U}*Conj[{I}] ] ; In DomainC_Mag ; }
          Term { [ {Us}*Conj[{Is}] ] ; In DomainS_Mag ; }
        }
      }

      { Name R ; Value {
          Term { [ -Re[{U}/{I}] ] ; In DomainC_Mag ; }
          Term { [ -Re[{Us}/{Is}] ] ; In DomainS_Mag ; }
        }
      }

      { Name L ; Value {
          Term { [ -Im[{U}/{I}]/(2*Pi*Freq) ] ; In DomainC_Mag ; }
          Term { [ -Im[{Us}/{Is}]/(2*Pi*Freq) ] ; In DomainS_Mag ; }
        }
      }

      // For convenience... possible scaling
      { Name R_per_km ; Value {
          Term { [ -Re[{U}/{I}]*1e3 ] ; In DomainC_Mag ; }
          Term { [ -Re[{Us}/{Is}]*1e3 ] ; In DomainS_Mag ; }
        }
      }

      { Name mL_per_km ; Value {
          Term { [ -1e6*Im[{U}/{I}]/(2*Pi*Freq) ] ; In DomainC_Mag ; }
          Term { [ -1e6*Im[{Us}/{Is}]/(2*Pi*Freq) ] ; In DomainS_Mag ; }
        }
      }

    }
  }

}


PostOperation{
  // Magnetic
  //-------------------------------

  { Name Mag_Maps; NameOfPostProcessing Magnetodynamic_a_2D;
    Operation {
      // You may want to see the result only in part of the domain => adapt domain that follows OnElementsOf
      // Name is not compulsory, it can be modified
      // Echo instructions are not necessary, useful for changing default visualisation options in the GUI

      Print[ az, OnElementsOf Domain_Mag, //Smoothing 1,
        Name "flux lines: Az [T m]", File "res/az.pos" ];
      Echo[Str["View[PostProcessing.NbViews-1].RangeType = 3;" ,// per timestep
          "View[PostProcessing.NbViews-1].NbIso = 25;" ,
          "View[PostProcessing.NbViews-1].IntervalsType = 1;" // isolines
        ], File "res/maps.opt"];

      Print[ b, OnElementsOf Domain_Mag, //Smoothing 1,
        Name "B [T]", File "res/b.pos" ];
      Echo[Str["View[PostProcessing.NbViews-1].RangeType = 3;" ,// per timestep
          "View[PostProcessing.NbViews-1].IntervalsType = 2;"
        ], File "res/maps.opt"];

      Print[ bm , OnElementsOf Domain_Mag,
        Name "|B| [T]", File "res/bm.pos" ];
       Echo[Str["View[PostProcessing.NbViews-1].RangeType = 3;" ,// per timestep
          "View[PostProcessing.NbViews-1].ShowTime = 0;",
          "View[PostProcessing.NbViews-1].IntervalsType = 2;"
        ], File "res/maps.opt"];

      Print[ jz , OnElementsOf Region[{DomainC_Mag, DomainS_Mag}],
        Name "jz [A/m^2] Conducting domain", File "res/jz_inds.pos" ];
      Echo[Str["View[PostProcessing.NbViews-1].RangeType = 3;" ,// per timestep
          "View[PostProcessing.NbViews-1].IntervalsType = 2;"
        ], File "res/maps.opt"];

      Print[ rhoj2 , OnElementsOf Region[{DomainC_Mag, DomainS_Mag}],
        Name "Power density", File "res/rhoj2.pos" ];
      Echo[Str["View[PostProcessing.NbViews-1].RangeType = 3;" ,// per timestep
          "View[PostProcessing.NbViews-1].ShowTime = 0;",
          "View[PostProcessing.NbViews-1].IntervalsType = 2;"
        ], File "res/maps.opt"];


      Print[ jm , OnElementsOf DomainC_Mag,
        Name "|j| [A/m^2] Conducting domain", File "res/jm.pos" ];
       Echo[Str["View[PostProcessing.NbViews-1].RangeType = 3;" ,// per timestep
          "View[PostProcessing.NbViews-1].ShowTime = 0;",
          "View[PostProcessing.NbViews-1].IntervalsType = 2;"
        ], File "res/maps.opt"];

    }
  }

  // This two definitions refer to the GUI => organisation of the results
  po = "{01Losses/";
  po2 = "{02PU-parameters/";
  { Name Mag_Global; NameOfPostProcessing Magnetodynamic_a_2D;
    Operation {
      // You may restrict DomainC_Mag to part of it

      Print[ JouleLosses[DomainC_Mag], OnGlobal, Format Table,
        SendToServer StrCat[po,"0Total conducting domain"],
        Units "W/m", File "res/losses_total.dat" ];

      Print[ JouleLosses[Inds], OnGlobal, Format Table,
        SendToServer StrCat[po,"3Source (stranded OR massive)"],
        Units "W/m", File "res/losses_inds.dat" ];

       Print[ JouleLosses[Steel], OnGlobal, Format Table,
        SendToServer StrCat[po,"1Steel"],
        Units "W/m", File "res/losses_steel.dat" ];

       Print[ JouleLosses[APLSheath], OnGlobal, Format Table,
        SendToServer StrCat[po,"2APLSheath"],
        Units "W/m", File "res/losses_APLSheath.dat" ];

      // Region to adapt according to your cable
      Print[ R, OnRegion Ind_1, Format Table,
        SendToServer StrCat[po2,"0R"],
        Units "Ω", File "res/Rinds.dat" ];
      Print[ L, OnRegion Ind_1, Format Table,
        SendToServer StrCat[po2,"1L"],
        Units "H", File "res/Linds.dat" ];

    }
  }
}
