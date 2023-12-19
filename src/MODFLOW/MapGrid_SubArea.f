C************************************************************************************************
C     Create linkage between MODFLOW Grid Cells and APEX subareas
C************************************************************************************************

      subroutine MapGrid_SubArea

      implicit none

C	Declare Variables
      integer   d,h,i,j,n,count,
     &          num_subarea,subarea_id,max_subarea,
     &          cell_id,num_lines,num_cells,line_count,dum,max_cells,
     &          nrow,ncol
      
      real      subarea_area,hru_area,cell_area,overlap_area

      integer,dimension(:), allocatable::subarea_ids,cell_active,
     &                                   cell_num_subareas,
     &                                   subarea_active,
     &                                   subarea_num_cells,
     &                                   cell_row_ids,cell_col_ids
      real,dimension(:), allocatable::area_fraction


      !open files
      open(3,file='MODFLOW/link_grid_sa')  !Ali
      open(12,file='MODFLOW/apexmf_grid2sa.txt')  !Ali


      print *, '     Grid Cells --> Subareas...'

      read(3,*) num_lines
      read(3,*) num_subarea
      read(3,*)
      read(3,*)
      read(3,*) !skip header information

      !allocate arrays
      allocate(area_fraction(200000))

      !First: Determine which subareas are present
      allocate(subarea_active(num_subarea))
      subarea_active = 0
      do n=1,num_lines
        read(3,*) cell_id,cell_area,subarea_id
        subarea_active(subarea_id) = 1  
      enddo      

      allocate(cell_row_ids(50000))
      allocate(cell_col_ids(50000))

      !Second: Determine the maximum number of grid cells per subarea
      rewind(3)
      read(3,*)
      read(3,*)
      read(3,*)
      read(3,*)
      read(3,*)
      allocate(subarea_num_cells(num_subarea))
      subarea_num_cells = 0
      max_cells = 0
      line_count = 0
      do n=1,num_subarea
        if(subarea_active(n).eq.1) then !only if cell is present
          subarea_id = n
          count = 1
          do while (subarea_id.eq.n) 
            if(line_count.lt.num_lines) then
              read(3,*) cell_id,cell_area,subarea_id
              count = count + 1
              line_count = line_count + 1
            else
              count = count -1
              goto 40
            endif
          enddo
          count = count - 2
 40       subarea_num_cells(n) = count !number of subareas that overlap the cell
          !determine if the number of subareas is the current maximum
          if(count.gt.max_cells) then
            max_cells = count
          endif
          backspace(3) !go back one line
          line_count = line_count - 1
        endif
      enddo !go to the next cell

      !Third: Calculate area fractions for each cell, and write out results
      rewind(3)
      read(3,*)
      read(3,*)
      read(3,*) nrow
      read(3,*) ncol
      read(3,*)
      write(12,101) num_subarea,max_cells
      do n=1,num_subarea 
        write(12,101) n,subarea_num_cells(n)
        if(subarea_active(n).eq.1) then
          !loop through the cells in the subarea
          do i=1,subarea_num_cells(n)
            read(3,*) cell_id,cell_area,subarea_id,
     &                overlap_area,subarea_area
            area_fraction(i) = overlap_area/subarea_area
            cell_row_ids(i) = int(cell_id/ncol) + 1
            cell_col_ids(i) = cell_id - (ncol * (cell_row_ids(i)-1))
          enddo
          !write out information for the subarea
          write(12,101) (cell_row_ids(i),i=1,subarea_num_cells(n))
          write(12,101) (cell_col_ids(i),i=1,subarea_num_cells(n))
          write(12,100) (area_fraction(i),i=1,subarea_num_cells(n))
        else
          write(12,*)
          write(12,*)
          write(12,*)
        endif
      enddo


 100  format(5000(f13.5))
 101  format(5000(i13))
    
      close(3)
      close(12)
      
	end
