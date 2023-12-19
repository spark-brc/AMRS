C************************************************************************************************
C     Create linkage between MODFLOW River Cells and APEX Subareas
C************************************************************************************************

      subroutine MapRiver_Grid

      implicit none

C	Declare Variables
      integer   d,h,i,j,k,n,count,
     &          cell_id,num_lines,num_cells,line_count,dum,
     &          num_river_cells,river_cell_1,river_cell_2,grid_id

      integer,dimension(:), allocatable::subbasin
      real,dimension(:), allocatable::river_length



      !open files
      open(4,file='MODFLOW/link_river_grid')  !Ali
      open(13,file='MODFLOW/apexmf_river2grid.txt')  !Ali

      print *, '     Subareas  <--> RiverCells...'


      !First: Determine the number of River Cells
      read(4,*) num_lines
      read(4,*) !skip header information

      num_river_cells = 0 
      do n=1,num_lines
        read(4,*) river_cell_1
        if(n.lt.num_lines) then
          read(4,*) river_cell_2
          if(river_cell_1.ne.river_cell_2) then
            num_river_cells = num_river_cells + 1
          endif
          backspace(4)
        else !for last line (always should be counted)
          num_river_cells = num_river_cells + 1
        endif
      enddo
      write(13,*) num_river_cells

      !Second: determine the subbasin for each River Cell, and any fractions (if more than one subbasin
      !intersects a grid cell)
      rewind(4)
      read(4,*)
      read(4,*)
      allocate(subbasin(10000))
      allocate(river_length(10000))
      line_count = 0
      do n=1,num_river_cells  
        read(4,*) grid_id !ID of current River Cell
        backspace(4)
        count = 1
        cell_id = grid_id
        !loop through the subbasins in the River Cell
        do while (cell_id.eq.grid_id)
          if(line_count.lt.num_lines) then
            read(4,*) cell_id,subbasin(count),river_length(count)
            count = count + 1
            line_count = line_count + 1
          else
            count = count - 1
            goto 10
          endif
        enddo
        count = count - 2
        !write out to file
 10     write(13,101) n,grid_id,count
        write(13,101) (subbasin(i),i=1,count)
        write(13,100) (river_length(i),i=1,count)
        subbasin = 0
        river_length = 0.
        backspace(4)
        line_count = line_count - 1
      enddo !go to next river cell

 
 100  format(5000(f13.5))
 101  format(5000(i13))
      
      close(4)   !Ali
      close(13)  !Ali
    
	end
