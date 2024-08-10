####### string manipulation #######

def add_spaces(string, length = 1, where = "left"):
    spaces = " " * (length - len(string))
    if where == "left":
        string = spaces + string
    elif where == "right":
        string += spaces
    return string
    
def str_round(num, dec):
    num_str = str(round(num, dec))
    l = len(num_str)
    li = len(str(int(num))) + dec + 1
    if l < li:
        num_str += "0" * (li - l)
    return num_str

def round_list(data, dec = 1):
    data = list(data)
    for i in range(len(data)):
        t = str(type(data[i]))
        if "float" in t or "int" in t:
            data[i] = str_round(float(data[i]), dec)
    return data

def list_to_string(data, where = "left", dec = 1):
    data = round_list(data, dec)
    data = list(map(str, data))
    spaces = max(map(len, data))
    return [add_spaces(el, spaces, where) for el in data]

colors = ['norm', 'black', 'gray', 'red', 'green', 'yellow', 'orange', 'blue', 'violet', 'cyan', 'bold']
color_codes = [0, 30, 2, 91, 92, 93, 33, 94, 95, 96, 1]
    
# it applies the proper color codes to a string
def set_color(text = "", color = "norm"):
    code = '\033['
    if type(color) == str:
        for c in range(len(colors)):
            if color == colors[c]:
                code += str(color_codes[c])
    code += 'm'
    return code + text + '\033[0m'

def print_snr(snr): # signal to noise print
    string = "fit snr: " + str(round(snr, 1))
    if snr < 10:
        print(set_color(string, "red"))
    elif snr < 30:
        print(set_color(string, "orange"))
    else:
        print(set_color(string, "green"))