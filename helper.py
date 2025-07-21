def sgn(x):
    if x == 0:
        return 0
    elif x < 0:
        return -1
    else:
        return 1

def map_range(x, input_start, input_end, output_start, output_end):
    return  (x - input_start) / (input_end - input_start) * (output_end - output_start) + output_start
