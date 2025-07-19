! ##################################################################################################################################
! Begin MIT license text.                                                                                    
! _______________________________________________________________________________________________________
                                                                                                         
! Copyright 2022 Dr William R Case, Jr (mystransolver@gmail.com)                                              
                                                                                                         
! Permission is hereby granted, free of charge, to any person obtaining a copy of this software and      
! associated documentation files (the "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to   
! the following conditions:                                                                              
                                                                                                         
! The above copyright notice and this permission notice shall be included in all copies or substantial   
! portions of the Software and documentation.                                                                              
                                                                                                         
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS                                
! OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,                            
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE                            
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER                                 
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,                          
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN                              
! THE SOFTWARE.                                                                                          
! _______________________________________________________________________________________________________
                                                                                                        
! End MIT license text.                                                                                      

      MODULE MITC_STUF

! This module contains variables that are calculated once for each element and used in various places.
  
      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE
  
      IMPLICIT NONE

      SAVE
    

! **********************************************************************************************************************************

                                                           ! Maximum number of nodes for an element
      INTEGER(LONG), PARAMETER        :: MELGP = 8         

                                                           ! Director vector at each element node
      REAL(DOUBLE)                    :: DIRECTOR(MELGP,3)  

                                                           ! Thickness in the direction of the director vector at element nodes.
                                                           ! Called a_k in ref [2] where k is node number.
      REAL(DOUBLE)                    :: DIR_THICKNESS(MELGP)

                                                           ! R and S coordinates of each element node
      REAL(DOUBLE)                    :: GP_RS(2,MELGP)


! **********************************************************************************************************************************

      END MODULE MITC_STUF
