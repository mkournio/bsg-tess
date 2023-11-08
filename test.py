

# declare a class 
class Employee:  
    
    # class attribute 
    count = 0      
        
    # define a method 
    def increase(self):  
        Employee.count += 1
    
# create an Employee  
# class object 
a1 = Employee()  
  
# calling object's method 
a1.increase()  
  
# print value of class attribute 
print(a1.count)  
     
a2 = Employee()  
  
a2.increase()  
  
print(a2.count)  
    
print(Employee.count)

