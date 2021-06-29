def kronecker_delta(a, b):
    '''
    Method to calculate the Kronecker delta.

            Parameters:
                    a (int): First index of Kronecker delta
                    b (int): Second index of Kronecker delta

            Returns:
                    (int): 0 or 1 as specified by the definition of the Kronecker delta
    '''
    
    if a == b:
        return 1
    else:
        return 0