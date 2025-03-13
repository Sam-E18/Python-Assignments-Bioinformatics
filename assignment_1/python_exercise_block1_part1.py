######Assignment_1:  There are four exercise to get in touch with python.
# Four initial exercises to begin to see the main components of Python. #######

#1)##Write a function that returns a float corresponding to the volume of a sphere##

def get_sphere_volume(radius):
    """ Return the volume of a sphere"""
    if radius < 0:
        return "Error value, cannot be negative! :)"
    pi = 3.1415  # I only use 4 decimals
    v_formula = (4 / 3) * pi * (radius ** 3)  
    return v_formula
print(get_sphere_volume(radius=5))

#Results: 523.5833333333333

##Other ways of performing the same exercise.##
import math 
##PD:Python modules are usually added at the beginning of the script. In this case I am adding it here to differentiate it from the previous one. ;) .
def get_sphere_volume(radius):
    """Return the volume of a sphere using the module math"""
    return (4/3) * math.pi * radius**3

get_sphere_volume(5)

#Results: 523.5987755982989

##2)Write a function that calculates and returns an integer corresponding to the factorial of an integer (n)###

#a)Using Recursivity

def recursive_factorial(n):
    """Return the factorial of the number"""
    if n < 0:
        return "Error: negative number"
    if n == 1:  
        return 1
    return n * recursive_factorial(n - 1)  

print(recursive_factorial(4))  

#Results: 24



#b)Without recursivity
def factorial(n):
    """Return the factorial of the number"""
    if n < 0:
        return "Error: negative number"
    number = 1
    for i in range(2, n + 1):  
        number *= i
    return number

print(factorial(5))  
#Results: 120

##Other ways of performing the same exercise.##
def factorial(n):
    
    if n < 0:
        return
    
    f = 1
    while n > 0:
        f = f * n
        n = n - 1

    return f

factorial(5)
#Results: 120

#3) Write a function for counting up numbers from 0 to n, showing the count up in the screen. 
# If parameter odd is set to True, prints only odd numbers#####

#a)Recursivity

def recursive_count_up(n, odd = True):
    """Count numbers from 0 to n using recursive function, if odd is true, return only odds numbers"""
    if n < 0:
        return 
    recursive_count_up(n - 1, odd)
    if odd == True:
        if n % 2 == 1:
            print(n)
    else:
        print(n)
recursive_count_up(5, odd=True)

#Results: If you set FALSE , count the numbers to n, if you set TRUE, counts only the odds.




#b)without recursivity
def count_up(n, odd=False):
    """Count numbers from 0 to n, if odd is true, return only odds numbers"""
    if odd == True:
        for i in range(1, n + 1 , 2):   
            print(i)
    else:
        for i in range(0, n + 1):
            print(i)
count_up(5, odd = True)


#4)###resolve the debug #Find and solve the bugs in the following function:

# def get_final_price(discount_percentage=10, price):
# 	"""Return the final price after applying the discount percentage"""
# 	return ((price + price) * percentage) / 100


def get_final_price( price, discount_percentage=10,):
   """Return the final price after applying the discount percentage. The
    discount formula is finalprice= price-(price*discount/100) """
   return (price - (price*discount_percentage/100)) 
print(get_final_price (price = 25.4))