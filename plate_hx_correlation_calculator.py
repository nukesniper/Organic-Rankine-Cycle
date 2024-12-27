def get_heat_transfer_coefficients(beta, Re):
    """
    Calculate the coefficients n and C_h based on beta and Reynolds number.

    Parameters:
        beta (float): Angle in degrees (e.g., 30, 45, 50, etc.)
        Re (float): Reynolds number

    Returns:
        tuple: (C_h, n) coefficients
    """
    if beta <= 30:
        if Re <= 10:
            return 0.718, 0.349
        else:
            return 0.348, 0.663
    elif beta == 45:
        if Re < 10:
            return 0.718, 0.349
        elif 10 <= Re <= 100:
            return 0.400, 0.598
        else:
            return 0.300, 0.663
    elif beta == 50:
        if Re < 20:
            return 0.630, 0.333
        elif 20 <= Re <= 300:
            return 0.291, 0.591
        else:
            return 0.130, 0.732
    elif beta == 60:
        if Re < 20:
            return 0.562, 0.326
        elif 20 <= Re <= 200:
            return 0.306, 0.529
        else:
            return 0.108, 0.703
    elif beta >= 65:
        if Re < 20:
            return 0.562, 0.326
        elif 20 <= Re <= 500:
            return 0.331, 0.503
        else:
            return 0.087, 0.718
    else:
        raise ValueError("Invalid beta value. Please use one of the specified angles in the table.")