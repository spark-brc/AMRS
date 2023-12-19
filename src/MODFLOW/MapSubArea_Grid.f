C************************************************************************************************
C     Create linakge between APEX subareas and the MODFLOW Grid Cells
C************************************************************************************************

      subroutine MapSubArea_Grid

      implicit none

      !declare Variables
      integer   d,h,i,j,n,count,
     &          num_subarea,subarea_id,max_subarea,
     &          cell_id,num_lines,num_cells,line_count,dum,max_cells
      
      real      subarea_area,hru_area,cell_area,overlap_area

      integer,dimension(:), allocatable::subarea_ids,cell_active,
     &                                   cell_num_subareas,
     &                                   subarea_active
      !real,dimension(:), allocatable::area_fraction
      real,dimension(:), allocatable::area_sub


      !open files
      open(2,file='MODFLOW/link_sa_grid')
      open(11,file='MODFLOW/apexmf_sa2grid.txt')

      print *, '     Subareas   --> Grid Cells...'

      read(2,*) num_lines
      read(2,*) num_cells
      read(2,*) !skip header information

      !allocate arrays
      !allocate(area_fraction(200000))
      allocate(area_sub(200000))
      allocate(subarea_ids(200000))

      !First: Determine which grid cells intersect the subareas
      allocate(cell_active(num_cells))
      cell_active = 0
      do n=1,num_lines
        read(2,*) cell_id
        cell_active(cell_id) = 1  
      enddo

      !Second: Determine the maximum number of subareas per grid cell
      rewind(2)
      read(2,*)
      read(2,*)
      read(2,*)
      allocate(cell_num_subareas(num_cells))
      cell_num_subareas = 0
      max_subarea = 0
      line_count = 0
      do n=1,num_cells
        if(cell_active(n).eq.1) then !only if cell is present
          cell_id = n
          count = 1
          do while (cell_id.eq.n) 
            if(line_count.lt.num_lines) then
              read(2,*) cell_id,cell_area,subarea_id
              count = count + 1
              line_count = line_count + 1
            else
              count = count -1
              goto 30
            endif
          enddo
          count = count - 2
 30       cell_num_subareas(n) = count !number of subareas that overlap the cell
          !determine if the number of subareas is the current maximum
          if(count.gt.max_subarea) then
            max_subarea = count
          endif
          backspace(2) !go back one line
          line_count = line_count - 1
        endif
      enddo !go to the next cell

      !Third: Calculate area fractions for each subarea, and write out results
      rewind(2)
      read(2,*)
      read(2,*)
      read(2,*)
      write(11,101) num_cells,max_subarea
      do n=1,num_cells 
        write(11,101) n,cell_num_subareas(n)
        if(cell_active(n).eq.1) then
          !loop through the subareas in the cell
          do i=1,cell_num_subareas(n)
            read(2,*) cell_id,cell_area,subarea_id,
     &                overlap_area,subarea_area
            area_sub(i) = overlap_area
            !area_fraction(i) = overlap_area/cell_area
            subarea_ids(i) = subarea_id  
          enddo
          !write out information for cell
          write(11,101) (subarea_ids(i),i=1,cell_num_subareas(n))
          !write(11,100) (area_fraction(i),i=1,cell_num_subareas(n))
          write(11,102) (area_sub(i),i=1,cell_num_subareas(n))
        endif
      enddo


 100  format(5000(f13.5))
 101  format(5000(i13))
 102  format(5000(f14.2))
      
      close(2)
      close(11)
    
	end
