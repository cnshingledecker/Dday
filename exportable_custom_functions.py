def is_int(val):
    try:
        int(val)
        return True
    except ValueError:
        return False

def is_float(val):
    try:
        float(val)
        return True
    except ValueError:
        return False

def find_nearest_index(value_to_find, pos, in_order_list):
    index = 0
    smallest_index = 0
    smallest_diff = 1e20
    for value in in_order_list:
        diff = abs(float(value_to_find) - float(value[pos]))
        if diff < smallest_diff:
            smallest_index = index
            smallest_diff = diff
        elif diff > smallest_diff:
            return smallest_index
        index += 1
    return smallest_index