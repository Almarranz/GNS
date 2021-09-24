FUNCTION RemoveRows, array, row

   ; array -- A 2D array from which rows will be removed.
   ; rows -- A vector of row indices to remove from array.

   ; Need both positional parameters.
   IF N_Params() NE 2 THEN BEGIN
       Print, "Usage: 'RemoveRows, array, rowsToRemove'"
       RETURN, -1
    ENDIF
    
    ; The array must be at least 2D.
    ndims = Size(array, /N_DIMENSIONS)
    IF ndims NE 2 THEN BEGIN
    
    stop
;        void = Dialog_Message('Array must be 2D.')
        Print, "Usage: 'RemoveRows, array, rowsToRemove'"
        RETURN, -1
    ENDIF
    
    ; The rows must be a vector.
    IF N_Elements(row) EQ 1 THEN rows = [row]
    
    ; Find the rows we are keeping.
    dims = Size(array, /DIMENSIONS)
    allrows = Lindgen(dims[1])
    goodRows = SetDifference(allrows, row)
    
    ; Return the subscripted array.
    RETURN, array[*,goodRows]

END