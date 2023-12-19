      subroutine salt_alloc
      !this subroutine reads data for RT3D's salinity package, and then allocates arrays
      
      use PARM, only:Ksp12,Ksp22,Ksp32,Ksp42,Ksp52,salt_solid_aqu
      use GLOBAL, only:NROW,NCOL,NLAY !MODFLOW
      use rt_global, only:insalt,saltdissolve_aqu,salt_aqu
      
      implicit none
      character (len=80) :: header
      integer i,j,k,m
      
      !read in Solubility Product parameters
      read(insalt,*) header
      read(insalt,*) Ksp12
      read(insalt,*) Ksp22
      read(insalt,*) Ksp32
      read(insalt,*) Ksp42
      read(insalt,*) Ksp52
      
      !read in aquifer salt mineral data (initial fractions)
      allocate(salt_solid_aqu(nrow,ncol,nlay,5))
      read(insalt,*) header
      do k=1,nlay
        !loop through the salt minerals
        read(insalt,*) header
        do m=1,5
          read(insalt,*) header  
          do i=1,nrow
            read(insalt,*) (salt_solid_aqu(i,j,k,m),j=1,ncol)
          enddo
        enddo
      enddo

      !other arrays
      allocate(saltdissolve_aqu(nrow,ncol,nlay))
      allocate(salt_aqu(nrow,ncol,nlay))
      
      return
      end

      
      
      
      subroutine salt_solve
      !this subroutine calculates salt ion concentrations based on equilibrium chemical reactions
	!(precipitation-dissolution, complexation, cation exchange)
             
      use parm
      use rt_global, only:cnew,prsity,saltdissolve_aqu,salt_aqu !RT3D
      use GLOBAL, only:NROW,NCOL,NLAY,HNEW,BOTM,IBOUND,DELC,DELR !MODFLOW
      
      implicit none

      integer i,j,k,m,n,dum,iter_count
      real    ion1,ion2,ion3,ion4,ion5,ion6,ion7,ion8,
     &        hru_area_m2,water_volume,gw_volume
      real    sol_water,sol_thick,waterC,
     &        Sol_CaCO3_p(1000),Sol_MgCO3_p(1000),Sol_CaSO4_p(1000),
     &        Sol_MgSO4_p(1000),Sol_NaCl_p(1000)
      real    I_Prep_in,I_diff,SkipedIEX
      real    cell_gw_volume,thickness,mass_before,mass_after,
     &        total_mass,mass_so4,mass_ca,mass_mg,mass_na,
     &				mass_k,mass_cl,mass_co3,mass_hco3,
     &        mass_before_dis,mass_before_sol,total_before,
     &        mass_after_dis,mass_after_sol,total_after,diff_total
      double precision IonStr,IS_temp,
     &                 K_ADJ1,K_ADJ2,K_ADJ3,K_ADJ4,K_ADJ5,
     &                 error1ST,error2ND,error3RD,errorTotal
      
      !track salt ion mass
      mass_before = 0.
      mass_after = 0.

      mass_before_dis = 0.
      mass_before_sol = 0.
      mass_after_dis = 0.
      mass_after_sol = 0.
      total_before = 0.
      total_after = 0.

      !loop through the grid cells - calculate change in concentration due to equilibrium reactions
      do k=1,nlay
        do i=1,nrow
          do j=1,ncol
            if(ibound(j,i,k).gt.0) then !only proceed if cell is active
              
              !volume of groundwater in the cell (m3)
              if(hnew(j,i,k).gt.botm(j,i,0)) then !head above the top of layer 1
                thickness = hnew(j,i,k) - botm(j,i,0)
              else
                if(hnew(j,i,k).lt.botm(j,i,k-1) .and. !BOTM(J,I,0) contains the top of layer 1
     &            hnew(j,i,k).gt.botm(j,i,k)) then 
                  thickness = hnew(j,i,k) - botm(j,i,k)
                else
                  thickness = botm(j,i,k-1) - botm(j,i,k)
                endif
              endif 
              cell_gw_volume = delr(j)*delc(i)*thickness*prsity(j,i,k)
                            
              !track salt ion mass in the cell
              mass_before = 0.
              mass_after = 0.  
              
              !check total salt mass (kg) in the cell
              mass_so4 = (cnew(j,i,k,3) * cell_gw_volume) / 1000.
              mass_ca = (cnew(j,i,k,4) * cell_gw_volume) / 1000.
              mass_mg = (cnew(j,i,k,5) * cell_gw_volume) / 1000.
              mass_na = (cnew(j,i,k,6) * cell_gw_volume) / 1000.
              mass_k = (cnew(j,i,k,7) * cell_gw_volume) / 1000.
              mass_cl = (cnew(j,i,k,8) * cell_gw_volume) / 1000.
              mass_co3 = (cnew(j,i,k,9) * cell_gw_volume) / 1000.
              mass_hco3 = (cnew(j,i,k,10) * cell_gw_volume) / 1000.
              mass_before = mass_so4 + mass_ca + mass_mg + mass_na +
     &										 mass_k + mass_cl + mass_co3 + mass_hco3
              
              !saturated water content (for saturated zone of aquifer) = porosity
              waterC = prsity(j,i,k)
              
              !retrieve the solid salt mineral concentration in the aquifer, and perform conversions
              Sol_CaCO3_p(1) = salt_solid_aqu(i,j,k,1)
              Sol_MgCO3_p(1) = salt_solid_aqu(i,j,k,2)
              Sol_CaSO4_p(1) = salt_solid_aqu(i,j,k,3)
              Sol_MgSO4_p(1) = salt_solid_aqu(i,j,k,4)
              Sol_NaCl_p(1)  = salt_solid_aqu(i,j,k,5)
              Sol_CaCO3(1) = (Sol_CaCO3_p(1)/100) *
     &                       (1.855/(waterC*100.0))*1000.0 
              Sol_MgCO3(1) = (Sol_MgCO3_p(1)/100) *
     &                       (1.855/(waterC*84.31))*1000.0 
              Sol_CaSO4(1) = (Sol_CaSO4_p(1)/100) *
     &                       (1.855/(waterC*136.14))*1000.0 
              Sol_MgSO4(1) = (Sol_MgSO4_p(1)/100) *
     &                       (1.855/(waterC*120.36))*1000.0 
              Sol_NaCl(1)  = (Sol_NaCl_p(1)/100) *
     &                       (1.855/(waterC*58.44))*1000.0 
              
              !retrieve the current (daily) salt ion solution concentrations from groundwater (mg/L) and convert to mol/L
              ion1 = cnew(j,i,k,3) !sulfate
              ion2 = cnew(j,i,k,4) !calcium
              ion3 = cnew(j,i,k,5) !magnesium
              ion4 = cnew(j,i,k,6) !sodium
              ion5 = cnew(j,i,k,7) !potassium
              ion6 = cnew(j,i,k,8) !chloride
              ion7 = cnew(j,i,k,9) !carbonate
              ion8 = cnew(j,i,k,10) !bicarbonate
              Sul_Conc(1) = ion1*((1.0/1000)*(1.0/96.06)) !sulfate
              Cal_Conc(1) = ion2*((1.0/1000)*(1.0/40.078)) !calcium
              Mg_Conc(1) = ion3*((1.0/1000)*(1.0/24.305)) !magnesium
              Sod_Conc(1) = ion4*((1.0/1000)*(1.0/23.0)) !sodium
              Pot_Conc(1) = ion5*((1.0/1000)*(1.0/39.0)) !potassium
              Cl_Conc(1) = ion6*((1.0/1000)*(1.0/35.45)) !chloride
              Car_Conc(1) = ion7*((1.0/1000)*(1.0/60.01)) !carbonate
              BiCar_Conc(1) = ion8*((1.0/1000)*(1.0/61.01)) !bicarbonate
     
              !define the activity coefficient using Extended Debye-Huckel equation
              call Ionic_strength(IS_temp,Cal_Conc(1),Sul_Conc(1),
     &                            Car_Conc(1),BiCar_Conc(1),Mg_Conc(1),
     &                            Sod_Conc(1),Pot_Conc(1))
              IonStr = IS_temp
              I_Prep_in = IonStr
              I_diff = 1
              
              !assign 1 for concentration to able count and store the data in while loop
              c11 = 1
              c22 = 1
              salt_c3 = 1
              salt_c4 = 1
              c5 = 1
              
              !This while loop compares I from precipitation loop with I from complexation. If they are about the same
              !(I_diff<0.001), it will go to the next row
              !do while (I_diff.ge.1e-2)
        
              call activity_coefficient(I_Prep_in)

              !update the K values 
              K_ADJ1 = LAMDA(1)*LAMDA(3)
              K_ADJ2 = LAMDA(5)*LAMDA(3)
              K_ADJ3 = LAMDA(1)*LAMDA(2)
              K_ADJ4 = LAMDA(5)*LAMDA(2)
              K_ADJ5 = LAMDA(6)*LAMDA(6) ! since I did not include the cl- yet

              salt_K1 = Ksp12/K_ADJ1
              salt_K2 = Ksp22/K_ADJ2
              salt_K3 = Ksp32/K_ADJ3
              salt_K4 = Ksp42/K_ADJ4
              salt_K5 = Ksp52/K_ADJ5

              errorTotal = 1

              !Precipitation-Dissolution package ----------------------------------------------------------------
              iter_count = 1
              do while (errorTotal.GE.1e-3)
                call CaCO3
                call MgCO3
                call CaSO4
                call MgSO4
                call NaCl

                !check the errors
                error1ST = Car_Conc(c22+1)-Car_Conc(c22+2)
                error2ND = Cal_Conc(c11+1)-Cal_Conc(c11+2)
                error3RD = Sul_Conc(salt_c4+1)-Sul_Conc(salt_c4+2)
                !errorTotal = ABS(MAX(error1ST,error2ND,error3RD))
             errorTotal = max(abs(error1ST),abs(error2ND),abs(error3RD))
        
                !update the counter for ions concentration
                c11 = c11 + 2
                c22 = c22 + 2
                salt_c3 = salt_c3 + 2
                salt_c4 = salt_c4 + 2
                c5 = c5 + 1

                !update iteration count
                iter_count = iter_count + 1
                if(iter_count.gt.500) then
                  goto 20
                endif

              enddo

              !convert mol/liter to ppm
 20           upion1 = Sul_Conc(salt_c4)
              upion2 = Cal_Conc(c11)
              upion3 = Mg_Conc(salt_c3)
              upion4 = Sod_Conc(c5)
              upion5 = Pot_Conc(1)
              upion6 = Cl_Conc(c5)
              upion7 = Car_Conc(c22)
              upion8 = BiCar_Conc(1)
        
              !enddo

              !update concentration after precipitation
              upion1 = Sul_Conc(salt_c4)*(96.06*1000.0)
              upion2 = Cal_Conc(c11)*(40.078*1000.0)
              upion3 = Mg_Conc(salt_c3)*(24.305*1000.0)
              upion4 = Sod_Conc(c5)*(23.0*1000.0)
              upion5 = Pot_Conc(1)*(39.0*1000.0)
              upion6 = Cl_Conc(c5)*(35.45*1000.0)
              upion7 = Car_Conc(c22)*(60.01*1000.0)
              upion8 = BiCar_Conc(1)*(61.01*1000.0) 

              !Cation Exchange package --------------------------------------------------------------------------
              call cationexchange

              !skipping cation exchange if: 
              if (upion2.le.0 .or. upion3.le.0 .or. upion4.le.0 .or. 
     &            upion5.le.0) then
                upion2 = Cal_Conc(c11)*(40.078*1000.0)
                upion3 = Mg_Conc(salt_c3)*(24.305*1000.0)
                upion4 = Sod_Conc(c5)*(23.0*1000.0)
                upion5 = Pot_Conc(1)*(39.0*1000.0)
                SkipedIEX = SkipedIEX + 1
              else
                upion2 = upion2
                upion3 = upion3
                upion4 = upion4
                upion5 = upion5
              endif

              !save data (upion = updated groundwater salt ion concentration)
              cnew(j,i,k,3) = upion1
              cnew(j,i,k,4) = upion2
              cnew(j,i,k,5) = upion3
              cnew(j,i,k,6) = upion4
              cnew(j,i,k,7) = upion5
              cnew(j,i,k,8) = upion6
              cnew(j,i,k,9) = upion7
              cnew(j,i,k,10) = upion8

              !check total salt mass (kg) in the cell
              mass_so4 = (cnew(j,i,k,3) * cell_gw_volume) / 1000.
              mass_ca = (cnew(j,i,k,4) * cell_gw_volume) / 1000.
              mass_mg = (cnew(j,i,k,5) * cell_gw_volume) / 1000.
              mass_na = (cnew(j,i,k,6) * cell_gw_volume) / 1000.
              mass_k = (cnew(j,i,k,7) * cell_gw_volume) / 1000.
              mass_cl = (cnew(j,i,k,8) * cell_gw_volume) / 1000.
              mass_co3 = (cnew(j,i,k,9) * cell_gw_volume) / 1000.
              mass_hco3 = (cnew(j,i,k,10) * cell_gw_volume) / 1000.
              mass_after = mass_so4 + mass_ca + mass_mg + mass_na +
     &										 mass_k + mass_cl + mass_co3 + mass_hco3

              !convert solids concentrations to solid percentage, and save
              Sol_CaCO3_p(c5) = (Sol_CaCO3(c5)*100.0)
     &                         /((1.855/(waterC*100.0))*1000.0) 
              Sol_MgCO3_p(c5) = (Sol_MgCO3(c5)*100.0)
     &                         /((1.855/(waterC*84.31))*1000.0) 
              Sol_CaSO4_p(c5) = (Sol_CaSO4(c5)*100.0)
     &                         /((1.855/(waterC*136.14))*1000.0) 
              Sol_MgSO4_p(c5) = (Sol_MgSO4(c5)*100.0)
     &                         /((1.855/(waterC*120.36))*1000.0) 
              Sol_NaCl_p(c5)  = (Sol_NaCl(c5)*100.0)
     &                         /((1.855/(waterC*58.44))*1000.0) 
              
              !save solids concentration for the HRU and the layer
              salt_solid_aqu(i,j,k,1) = Sol_CaCO3_p(c5)
              salt_solid_aqu(i,j,k,2) = Sol_MgCO3_p(c5)
              salt_solid_aqu(i,j,k,3) = Sol_CaSO4_p(c5)
              salt_solid_aqu(i,j,k,4) = Sol_MgSO4_p(c5)
              salt_solid_aqu(i,j,k,5) = Sol_NaCl_p(c5)      

              !store cell value for aquifer salt mass balance (solid --> dissolved)
              saltdissolve_aqu(i,j,k) = mass_after - mass_before
              salt_aqu(i,j,k) = mass_before
              
            endif
          enddo
        enddo
      enddo

      
      return
      end


      

      ! Calculate Ionic Strength of Water *******************************************************************
      subroutine Ionic_Strength(IS_temp,A,B,C,D,E,F,G)
      
      implicit none

      double precision IS_temp,A,B,C,D,E,F,G
      real CharBal(7)
      
      DATA CharBal/2.0, -2.0, -2.0, -1.0, 2.0, 1.0, 1.0/     
      
      IS_temp = 0.5*(CharBal(1)**2*A
     &+CharBal(2)**2*B
     &+CharBal(3)**2*C
     &+CharBal(4)**2*D
     &+CharBal(5)**2*E
     &+CharBal(6)**2*F
     &+CharBal(7)**2*G)
     
      return
      end
      
      

      
      
      ! Calculate Activity Coefficient **********************************************************************
      subroutine activity_coefficient(I_Prep_in)
      
      use parm

      real CharBal(7),a_size(7),I_Prep_in
      DATA CharBal/2.0, -2.0, -2.0, -1.0, 2.0, 1.0, 1.0/
      DATA a_size/6.0, 4.0, 4.5, 4.0, 8.0, 4.5, 3.0/ 
      A = 0.5 !at 298 K
      B = 0.33 !at 298 K
      
      if (I_Prep_in.LE.1e-1) then
        do ii = 1,7
          LAMDA(ii)= 10.0**(-A*CharBal(ii)**2.0
     &    *(I_Prep_in**0.5/(1+B*a_size(ii)*I_Prep_in**0.5)))
        enddo
      elseif (I_Prep_in.GE.5) then
        I_Prep_in = 0.5
        do ii = 1,7
          LAMDA(ii)= 10.0**(-A*CharBal(ii)**2.0
     &    *(I_Prep_in**0.5/(1+I_Prep_in**0.5)-0.3*I_Prep_in))
        enddo
      else
        do ii = 1,7
          LAMDA(ii)= 10.0**(-A*CharBal(ii)**2.0
     &    *(I_Prep_in**0.5/(1+I_Prep_in**0.5)-0.3*I_Prep_in))
        enddo
      endif
      
      return
      end
           
      



      ! CaSO4 ***********************************************************************************************
      !disp('**************************************************************')
      !disp('The reaction is: CaSO4(s)<--> Ca2+(aq) + SO42-(aq)')
      !disp('The Ksp(Solubility Product Constants) of CaSO4(s) is 4.93*10^-5')
      !disp('**************************************************************') 
      subroutine CaSO4

      use parm

      implicit none

      Double Precision M1,M2,M3,Ksp,Solv,Trial_Ksp,
     & PosSolv,CalSul_Prep,Solid_CaSO4,Dissolved_Solid,
     & Calcium_Conc,Sulfate_Conc
      
      M1 = Sol_CaSO4(c5)
      M2 = Cal_Conc(c11+1)
      M3 = Sul_Conc(salt_c4)
      Ksp = salt_K3
      
      Solv = 0.5*(-(M2+M3)+sqrt((M2+M3)**2-4*(M2*M3-Ksp)))
      Trial_Ksp = M2 * M3

      if (Trial_Ksp.gt.Ksp) then
        !disp('precipitation WILL form')
        ! Defining temporary parameter PosSolve
        PosSolv = abs(Solv)
        CalSul_Prep = PosSolv 
        Solid_CaSO4 = M1+CalSul_Prep
        Dissolved_Solid = 0
        Calcium_Conc = M2-CalSul_Prep
        Sulfate_Conc = M3-CalSul_Prep
      
      elseif (M1.gt.Solv) then
        !disp('Dissolution WILL form')
        !Defining temporary parameter PosSolv
        PosSolv = Solv
      
        !disp('****************************************')
        !disp('Portion of the solid will be dissolved')
        !disp('****************************************')
        Calcium_Conc = (M2+PosSolv)
        Sulfate_Conc = (M3+PosSolv)
        Dissolved_Solid = PosSolv
        Solid_CaSO4 = (M1-PosSolv)
      
      else 
        !disp('****************************************')
        !disp('Solid will be compeletly dissolved')
        !disp('****************************************')
        Calcium_Conc = (M1+M2)
        Sulfate_Conc = (M1+M3)
        Dissolved_Solid = M1
        Solid_CaSO4 = 0
      endif
      
      Sol_CaSO4(c5+1) = Solid_CaSO4
      Cal_Conc(c11+2) = Calcium_Conc
      Sul_Conc(salt_c4+1) = Sulfate_Conc

      return
      end





      ! MgCO3 ***********************************************************************************************
      !disp('**************************************************************')
      !disp('The reaction is: MgCO3(s)<--> Mg2+(aq) + CO32-(aq)')
      !disp('The Ksp(Solubility Product Constants) of MgCO3(s) is 4.7937*10^-6')
      !disp('**************************************************************') 
      subroutine MgCO3

      use parm

      implicit none

      integer j

      Double Precision M1,M2,M3,Ksp,Solv,Trial_Ksp,
     & PosSolv,MgCar_Prep,Solid_MgCO3,Dissolved_Solid,
     & Mag_Conc,Carbonate_Conc
      
      M1 = Sol_MgCO3(c5)
      M2 = Mg_Conc(salt_c3)
      M3 = Car_Conc(c22+1)
      Ksp = salt_K2
      Solv = 0.5*(-(M2+M3)+sqrt((M2+M3)**2-4*(M2*M3-Ksp)))
      Trial_Ksp=M2*M3

      if (Trial_Ksp.GT.Ksp) then   
        !disp('precipitation WILL form')
        ! Defining temporary parameter PosSolve
        PosSolv = abs(Solv)
        MgCar_Prep = PosSolv 
        Solid_MgCO3 = M1+MgCar_Prep
        Dissolved_Solid = 0
        Mag_Conc = M2-MgCar_Prep
        Carbonate_Conc = M3-MgCar_Prep
      
      elseif(M1.GT.Solv) then   
        !disp('Dissolution WILL form')
        !Defining temporary parameter PosSolv
        PosSolv = Solv
        ! IF (M1.GT.PosSolv) THEN 
        !disp('****************************************')
        !disp('Portion of the solid will be dissolved')
        !disp('****************************************')
        Mag_Conc = (M2+PosSolv)
        Carbonate_Conc = (M3+PosSolv)
        Dissolved_Solid = PosSolv
        Solid_MgCO3 = (M1-PosSolv)
      
      else 
        !disp('****************************************')
        !disp('Solid will be compeletly dissolved')
        !disp('****************************************')
        Mag_Conc = (M1+M2)
        Carbonate_Conc = (M1+M3)
        Dissolved_Solid = M1
        Solid_MgCO3= 0
      
      endif

      Sol_MgCO3(c5+1) = Solid_MgCO3
      Mg_Conc(salt_c3+1) =   Mag_Conc
      Car_Conc(c22+2) =  Carbonate_Conc  
         
      return
      end





      ! NaCl ************************************************************************************************
      !disp('**************************************************************')
      !disp('The reaction is: NaCl(s)<--> Na+(aq) + Cl-(aq)')
      !disp('The Ksp(Solubility Product Constants) of NaCl(s) is 37.3')
      !disp('**************************************************************')  
      subroutine NaCl

      use parm

      implicit none

      Double Precision M1,M2,M3,Ksp,Solv,Trial_Ksp,
     & PosSolv,SodiumChloride_Prep,Solid_NaCl,Dissolved_Solid,
     & Sodium_Conc,Chloride_Conc
      
      M1 = Sol_NaCl(c5)
      M2 = Sod_Conc(c5)
      M3 = Cl_Conc(c5)
      Ksp = salt_K5
      
      Solv = 0.5*(-(M2+M3)+sqrt((M2+M3)**2-4*(M2*M3-Ksp)))
      Trial_Ksp=M2*M3

      if (Trial_Ksp.GT.Ksp) then    
        !disp('precipitation WILL form')
        ! Defining temporary parameter PosSolve
        PosSolv = abs(Solv)
        SodiumChloride_Prep = PosSolv 
        Solid_NaCl = M1+SodiumChloride_Prep
        Dissolved_Solid = 0
        Sodium_Conc = M2-SodiumChloride_Prep
        Chloride_Conc  = M3-SodiumChloride_Prep
      
      elseif(M1.GT.Solv) then
        !disp('Dissolution WILL form')
        !Defining temporary parameter PosSolv
        PosSolv = Solv
        !  IF (M1.GT.PosSolv) THEN 
        !disp('****************************************')
        !disp('Portion of the solid will be dissolved')
        !disp('****************************************')
        Sodium_Conc = (M2+PoSSolv)
        Chloride_Conc  = (M3+PosSolv)
        Dissolved_Solid = PosSolv
        Solid_NaCl = (M1-PosSolv)
      
      else
        !disp('****************************************')
        !disp('Solid will be compeletly dissolved')
        !disp('****************************************')
        Sodium_Conc = (M1+M2)
        Chloride_Conc  = (M1+M3)
        Dissolved_Solid = M1
        Solid_NaCl= 0
      
      endif

      Sol_NaCl(c5+1) = Solid_NaCl
      Sod_Conc(c5+1) = Sodium_Conc
      Cl_Conc(c5+1) = Chloride_Conc
      
      return
      end





      ! MgSO4 ***********************************************************************************************
      !disp('**************************************************************')
      !disp('The reaction is: MgSO4(s)<--> Mg2+(aq) + SO42-(aq)')
      !disp('The Ksp(Solubility Product Constants) of MgSO4(s) is ?')
       !disp('**************************************************************')    
      subroutine MgSO4

      use parm

      implicit none

      Double Precision M1,M2,M3,Ksp,Solv,Trial_Ksp,
     & PosSolv,MgSul_Prep,Solid_MgSO4,Dissolved_Solid,
     & Mag_Conc,Sulfate_Conc
      
      M1 = Sol_MgSO4(c5)
      M2 = Mg_Conc(salt_c3+1)
      M3 = Sul_Conc(salt_c4+1)
      Ksp = salt_K4
      
      Solv = 0.5*(-(M2+M3)+sqrt((M2+M3)**2-4*(M2*M3-Ksp)))
      Trial_Ksp=M2*M3

      if (Trial_Ksp.GT.Ksp) then  
        !disp('precipitation WILL form')
        ! Defining temporary parameter PosSolve
        PosSolv = abs(Solv)
        MgSul_Prep = PosSolv 
        Solid_MgSO4 = M1+MgSul_Prep
        Dissolved_Solid = 0
        Mag_Conc = M2-MgSul_Prep
        Sulfate_Conc = M3-MgSul_Prep

      elseif(M1.GT.Solv) then
        !disp('Dissolution WILL form')
        !Defining temporary parameter PosSolv
        PosSolv = Solv
        !IF (M1.GT.PosSolv) THEN 
        !disp('****************************************')
        !disp('Portion of the solid will be dissolved')
        !disp('****************************************')
        Mag_Conc = (M2+PosSolv)
        Sulfate_Conc = (M3+PosSolv)
        Dissolved_Solid = PosSolv
        Solid_MgSO4 = (M1-PosSolv)

      else
        !disp('****************************************')
        !disp('Solid will be compeletly dissolved')
        !disp('****************************************')
        Mag_Conc = (M1+M2)
        Sulfate_Conc = (M1+M3)
        Dissolved_Solid = M1
        Solid_MgSO4= 0
      
      endif
        
      Sol_MgSO4(c5+1) = Solid_MgSO4
      Mg_Conc(salt_c3+2) = Mag_Conc
      Sul_Conc(salt_c4+2) = Sulfate_Conc
      
      return
      end





      ! CaCO3 ***********************************************************************************************
      !disp('**************************************************************')
      !disp('The reaction is: CaCO3(s)<--> Ca2+(aq) + CO32-(aq)')
      !disp('The Ksp(Solubility Product Constants) of CaCO3(s) is 3.0702*10^-9')
      !disp('**************************************************************')     
      subroutine CaCO3
      
      use parm
       
      Double Precision M1,M2,M3,Ksp,Solv,Trial_Ksp,
     & PosSolv,CalCar_Prep,Solid_CaCO3,Dissolved_Solid,
     & Calcium_Conc,Carbonate_Conc
            
      M1 = Sol_CaCO3(c5)
      M2 = Cal_Conc(c11)
      M3 = Car_Conc(c22)
      Ksp = salt_K1
      Solv = 0.5*(-(M2+M3)+sqrt((M2+M3)**2-4*(M2*M3-Ksp)))
      Trial_Ksp=M2*M3

      if (Trial_Ksp.GT.Ksp) then
        !disp('precipitation WILL form')
        ! Defining temporary parameter PosSolve
        PosSolv = abs(Solv)
        CalCar_Prep = PosSolv 
        Solid_CaCO3 = M1+CalCar_Prep
        Dissolved_Solid = 0
        Calcium_Conc = M2-CalCar_Prep
        Carbonate_Conc = M3-CalCar_Prep

      elseif(M1.GT.Solv) then
        !disp('Dissolution WILL form')
        !Defining temporary parameter PosSolv
        PosSolv = Solv
        !IF (M1.GT.PosSolv) THEN 
        !disp('****************************************')
        !disp('Portion of the solid will be dissolved')
        !disp('****************************************')
        Calcium_Conc = (M2+Solv)
        Carbonate_Conc = (M3+Solv)
        Dissolved_Solid = Solv
        Solid_CaCO3 = (M1-Solv)

      else
        !disp('****************************************')
        !disp('Solid will be compeletly dissolved')
        !disp('****************************************')
        Calcium_Conc = (M1+M2)
        Carbonate_Conc = (M1+M3)
        Dissolved_Solid = M1
        Solid_CaCO3= 0

      endif

      Sol_CaCO3(c5+1) = Solid_CaCO3 
      Cal_Conc(c11+1) =  Calcium_Conc
      Car_Conc(c22+1) = Carbonate_Conc 
      
      return
      end





      ! Calculate Cation Exchange ***************************************************************************   
      ! Developed by Saman Tavakoli 
      ! At Colorado State University 
      ! October 2016
      ! Includes Ca2+, Mg2+, Na+,K+
      ! CEC = Cation Exchange Capacity in meq/100g soil
      !Assumption for CEC: constant for a given soil, independent of of pH, ion type and concentration
      !Sel_K1 through Sel_K6 stands for selectivity coefficient of exchange reaction
      !XCAINI and others stands for inital ion which attached to the soil particle
      subroutine cationexchange
      
      use parm, only: upion2,upion3,upion4,upion5

      implicit none 

      real CEC,Sel_K1,Sel_K2,Sel_K3,Sel_K4,Sel_K5,Sel_K6,
     &     XCAINI,XMGINI,XNAINI,XKINI,
     &     DeltaX_Ca,DeltaX_Mg,DeltaX_Na,DeltaX_K,
     &     Con_Ca,Con_Mg,Con_Na,Con_K,
     &     X_Ca,X_Mg,X_Na,X_K
      
      !CEC selected based on soil type; for simplicity, for now used one value based on the sandy-loam soil type
      CEC  = 15 ! meq/100g soil 
      
      !Convert the ions concentration in solution from ppm to mmol/liter 
      Con_Ca = upion2/40  ! mmol/liter water 
      Con_Mg = upion3/24  ! mmol/liter water 
      Con_Na = upion4/23  ! mmol/liter water 
      Con_K  = upion5/39  ! mmol/liter water 

      !using Gapon Convention 
      Sel_K1 = 0.7
      Sel_K2 = 6
      Sel_K3 = 0.4
      Sel_K4 = 0.2
      Sel_K5 = 4
      Sel_K6 = 16

      !Make sure an appropriate value to include
      XCAINI = 13.6  ! meg/100g soil
      XMGINI = 0.17  !meq/100g soil
      XNAINI = 0.25 !meq/100g soil
      XKINI = 0.18 !meq/100g soil
      X_Ca = CEC/(1+((Sel_K1*(Con_Mg)**0.5)/((Con_Ca)**0.5))+
     &       ((Con_Na)/(Sel_K2*(Con_Ca)**0.5))+
     &       ((Con_K)/(Sel_K3*(Con_Ca)**0.5)))
      X_Mg = CEC/(((Con_Ca)**0.5/(Sel_K1*(Con_Mg)**0.5))+1+
     &       ((Con_Na)/(Sel_K4*(Con_Mg)**0.5))+
     &       ((Con_K)/(Sel_K5*(Con_Mg)**0.5)))
      X_Na = CEC/(((Sel_K2*(Con_Ca)**0.5)/(Con_Na))+
     &((Sel_K5*(Con_Mg)**0.5)/(Con_Na))+1+((Sel_K6*Con_K)/(Con_Na)))
      X_K = CEC/(((Sel_K3*(Con_Ca)**0.5)/(Con_K))+
     &((Sel_K4*(Con_Mg)**0.5)/(Con_K))+(Con_Na/(Sel_K6*Con_K))+1)
      
      !Post processing the solution concentration
      
C     Ref :http://www.public.iastate.edu/~teloynac/354ppcecsol.html 
C     1 meq Ca = 20 mg 
C     1 meq Mg = 12 mg 
C     1 meq Na = 23 mg 
C     1 meq K = 39 mg
      
C ppm = mg/kg 
C mg/ 100g soil * 10 = mg/kg 
C (X meq Ca/100g soil) * (20 mg/mequ)* 10 = 200X mg/kg or ppm
C (X meq Mg/100g soil) * (12 mg/mequ)* 10 = 120X mg/kg or ppm
C (X meq Na/100g soil) * (23 mg/mequ)* 10 = 230X mg/kg or ppm
C (X meq K/100g soil) * (39 mg/mequ)* 10 = 390X mg/kg or ppm
      
      
C Need Soil Bulk Density, assume : 1650 kg/m3
C Need water content from RT3D, for simplicity right now, water content = 0.3
      
      DeltaX_Ca = X_Ca - XCAINI
      DeltaX_Mg = X_Mg - XMGINI
      DeltaX_Na = X_Na - XNAINI
      DeltaX_K = X_K - XKINI
      
      upion2 = Con_Ca*40 - DeltaX_Ca*200*1.65/0.3 ! 0.3 = water content 
      upion3 = Con_Mg*24 - DeltaX_Mg*120*1.65/0.3
      upion4 = Con_Na*23 - DeltaX_Na*230*1.65/0.3
      upion5 = Con_K*39 - DeltaX_K*390*1.65/0.3

      if (Con_Ca.LE.0) then
          upion2 = -10
      endif 
      if (Con_Mg.LE.0) then
          upion3 = -10
      endif
      if (Con_Na.LE.0) then
          upion4 = -10
      endif
      if (Con_K.LE.0) then
          upion5 = -10
      endif
      
      return 
      end ! end subroutine cation exchange 
      