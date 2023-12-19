      subroutine salt_chem_soil !rtb
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine calculates salt ion concentrations based on equilibrium chemical reactions
!!    (precipitation-dissolution, complexation, cation exchange)
             
      use parm

      implicit none

      integer j,jj,m,n,dum,iter_count,LD1
      real    ion1,ion2,ion3,ion4,ion5,ion6,ion7,ion8,
     &        sub_area_m2,water_volume,gw_volume
      real    sol_water,sol_thick,waterC,
     &        Sol_CaCO3_p(1000),Sol_MgCO3_p(1000),Sol_CaSO4_p(1000),
     &        Sol_MgSO4_p(1000),Sol_NaCl_p(1000)
      real    I_Prep_in,I_diff,SkipedIEX
      real    soil_volume,salt_fract,mass_before(8),mass_after(8),
     &        salt_mass_kg,soil_mass,
     &        mass_before_dis,mass_before_sol,total_before,
     &        mass_after_dis,mass_after_sol,total_after,diff_total
      double precision IonStr,IS_temp,
     &                 K_ADJ1,K_ADJ2,K_ADJ3,K_ADJ4,K_ADJ5,
     &                 error1ST,error2ND,error3RD,errorTotal
      

      !subarea ID = ISA

      !area of the subarea in m2
      sub_area_m2 = WSA(ISA) * 10000.

      !track salt ion mass before and after equilibrium reactions
      mass_before(8) = 0.
      mass_after(8) = 0.
      mass_before_dis = 0.
      mass_before_sol = 0.
      mass_after_dis = 0.
      mass_after_sol = 0.
      total_before = 0.
      total_after = 0.
      
      !loop through the soil layers in the current subarea
      do jj=1,NBSL(ISA)
        
        !soil layer ID
        ISL = LID(jj,ISA)
        
        !calculate water content of the soil layer
        sol_thick = Z(ISL,ISA) !thickness of soil layer (mm)
        sol_water = SWST(ISL,ISA) !water in soil layer (mm)
        waterC = sol_water / sol_thick
        
        !retrieve the solid salt mineral concentration
        Sol_CaCO3_p(1) = salt_solid_soil(ISL,ISA,1)
        Sol_MgCO3_p(1) = salt_solid_soil(ISL,ISA,2)
        Sol_CaSO4_p(1) = salt_solid_soil(ISL,ISA,3)
        Sol_MgSO4_p(1) = salt_solid_soil(ISL,ISA,4)
        Sol_NaCl_p(1)  = salt_solid_soil(ISL,ISA,5)
        
        !calculate total salt mass from salt minerals
        soil_volume = sub_area_m2 * (sol_thick/1000.) !m3 of soil
        soil_mass = soil_volume * (BD(ISL,ISA)*1000.) !kg of soil
        do m=1,5
          mass_before_sol = mass_before_sol + 
     &          (soil_mass*(salt_solid_soil(ISL,ISA,m)/100.))
        enddo

        !perform conversions
        Sol_CaCO3(1) = (Sol_CaCO3_p(1)/100) *
     &                 (BD(ISL,ISA)/(waterC*100.0))*1000.0 
        Sol_MgCO3(1) = (Sol_MgCO3_p(1)/100) *
     &                 (BD(ISL,ISA)/(waterC*84.31))*1000.0 
        Sol_CaSO4(1) = (Sol_CaSO4_p(1)/100) *
     &                 (BD(ISL,ISA)/(waterC*136.14))*1000.0 
        Sol_MgSO4(1) = (Sol_MgSO4_p(1)/100) *
     &                 (BD(ISL,ISA)/(waterC*120.36))*1000.0 
        Sol_NaCl(1)  = (Sol_NaCl_p(1)/100) *
     &                 (BD(ISL,ISA)/(waterC*58.44))*1000.0 
        
        !get concentration of salt ions in the layer
        water_volume = (SWST(ISL,ISA)/1000.) * sub_area_m2 !m * m2 = m3
        do m=1,mion
          salt_mass_kg = wsalt(ISL,ISA,m) * WSA(ISA) !kg/ha * ha = kg of salt
          mass_before(m) = mass_before(m) + salt_mass_kg
          soil_csalt(ISL,ISA,m) = (salt_mass_kg * 1000.) / water_volume !g/m3 = mg/L
        enddo

        !total dissolved mass
        do m=1,mion
          mass_before_dis = mass_before_dis + mass_before(m)
        enddo

        !total salt mass in soil layer
        total_before = total_before + (mass_before_dis+mass_before_sol)

        !retrieve the current (daily) salt ion solution concentrations from soil water (mg/L) and convert to mol/L
        ion1 = soil_csalt(ISL,ISA,1) !sulfate
        ion2 = soil_csalt(ISL,ISA,2) !calcium
        ion3 = soil_csalt(ISL,ISA,3) !magnesium
        ion4 = soil_csalt(ISL,ISA,4) !sodium
        ion5 = soil_csalt(ISL,ISA,5) !potassium
        ion6 = soil_csalt(ISL,ISA,6) !chloride
        ion7 = soil_csalt(ISL,ISA,7) !carbonate
        ion8 = soil_csalt(ISL,ISA,8) !bicarbonate
        Sul_Conc(1) = ion1*((1.0/1000)*(1.0/96.06)) !sulfate
        Cal_Conc(1) = ion2*((1.0/1000)*(1.0/40.078)) !calcium
        Mg_Conc(1) = ion3*((1.0/1000)*(1.0/24.305)) !magnesium
        Sod_Conc(1) = ion4*((1.0/1000)*(1.0/23.0)) !sodium
        Pot_Conc(1) = ion5*((1.0/1000)*(1.0/39.0)) !potassium
        Cl_Conc(1) = ion6*((1.0/1000)*(1.0/35.45)) !chloride
        Car_Conc(1) = ion7*((1.0/1000)*(1.0/60.01)) !carbonate
        BiCar_Conc(1) = ion8*((1.0/1000)*(1.0/61.01)) !bicarbonate

        !define the activity coefficient using Extended Debye-Huckel equation
        call Ionic_strength(IS_temp,Cal_Conc(1),Sul_Conc(1),Car_Conc(1),
     &                      BiCar_Conc(1),Mg_Conc(1),Sod_Conc(1),
     &                      Pot_Conc(1))
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

          salt_K1 = Ksp11/K_ADJ1
          salt_K2 = Ksp21/K_ADJ2
          salt_K3 = Ksp31/K_ADJ3
          salt_K4 = Ksp41/K_ADJ4
          salt_K5 = Ksp51/K_ADJ5

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
              goto 10
            endif

          enddo 
          !************************ End Precipitation-Dissolution Package *************************

          !convert mol/liter to ppm
10        upion1 = Sul_Conc(salt_c4)
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
        upion8= BiCar_Conc(1)*(61.01*1000.0) 

        !Cation Exchange package --------------------------------------------------------------------------
        call cationexchange

        !skipping cation exchange if: 
        if (upion2.le.0 .or. upion3.le.0 .or. upion4.le.0 .or. 
     &      upion5.le.0) then
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

        !save data (upion = updated soil water salt ion concentration)
        soil_csalt(ISL,ISA,1) = upion1
        soil_csalt(ISL,ISA,2) = upion2
        soil_csalt(ISL,ISA,3) = upion3
        soil_csalt(ISL,ISA,4) = upion4
        soil_csalt(ISL,ISA,5) = upion5
        soil_csalt(ISL,ISA,6) = upion6
        soil_csalt(ISL,ISA,7) = upion7
        soil_csalt(ISL,ISA,8) = upion8

        !convert to kg/ha
        do m=1,mion
          wsalt(ISL,ISA,m) = (soil_csalt(ISL,ISA,m)/1000.)*water_volume
     &                       / WSA(ISA)    
        enddo

        !check mass
        do m=1,mion !total salt mass in soil layer
          salt_mass_kg = wsalt(ISL,ISA,m) * WSA(ISA) !kg of salt
          mass_after(m) = mass_after(m) + salt_mass_kg
        enddo

        !total dissolved mass
        do m=1,mion
          mass_after_dis = mass_after_dis + mass_after(m)
        enddo

        !convert solids concentrations to solid percentage, and save
        Sol_CaCO3_p(c5) = (Sol_CaCO3(c5)*100.0)
     &                    /((BD(ISL,ISA)/(waterC*100.0))*1000.0) 
        Sol_MgCO3_p(c5) = (Sol_MgCO3(c5)*100.0)
     &                    /((BD(ISL,ISA)/(waterC*84.31))*1000.0) 
        Sol_CaSO4_p(c5) = (Sol_CaSO4(c5)*100.0)
     &                    /((BD(ISL,ISA)/(waterC*136.14))*1000.0) 
        Sol_MgSO4_p(c5) = (Sol_MgSO4(c5)*100.0)
     &                    /((BD(ISL,ISA)/(waterC*120.36))*1000.0) 
        Sol_NaCl_p(c5)  = (Sol_NaCl(c5)*100.0)
     &                    /((BD(ISL,ISA)/(waterC*58.44))*1000.0) 
           
        !save solids concentration for the HRU and the layer
        salt_solid_soil(ISL,ISA,1) = Sol_CaCO3_p(c5)
        salt_solid_soil(ISL,ISA,2) = Sol_MgCO3_p(c5)
        salt_solid_soil(ISL,ISA,3) = Sol_CaSO4_p(c5)
        salt_solid_soil(ISL,ISA,4) = Sol_MgSO4_p(c5)
        salt_solid_soil(ISL,ISA,5) = Sol_NaCl_p(c5)
        
        !calculate total salt mass from salt minerals
        do m=1,5
          mass_after_sol = mass_after_sol + 
     &          (soil_mass*(salt_solid_soil(ISL,ISA,m)/100.))
        enddo

        !total salt mass in soil layer
        total_after = total_after + (mass_after_dis + mass_after_sol)

      enddo !go to next soil layer

      !store HRU value for soil salt (solid --> dissolved)
      do m=1,mion
        mass_before(m) = mass_before(m) / WSA(ISA) !kg/ha
        mass_after(m) = mass_after(m) / WSA(ISA) !kg/ha
        saltdissolve_soil(ISA,m) = mass_after(m) - mass_before(m)
      enddo
      
      !for the first soil layer, always set to initial salt mineral fraction
      !(in recognition that likely there is much more salt available then what is simulated; i.e. there is limitless salt reservoir in the top layer)
      LD1 = LID(1,ISA)
      do m=1,5
        salt_solid_soil(LD1,ISA,m) = init_salt_fraction(LD1,ISA,m)
      enddo
      
        
      return
      end

