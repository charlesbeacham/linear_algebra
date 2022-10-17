from math import acos, pi, sqrt
from decimal import Decimal, getcontext

getcontext().prec = 30

class Vector(object):
    def __init__(self, coordinates):
        try:
            if not coordinates:
                raise ValueError
            self.coordinates = tuple([Decimal(x) for x in coordinates])
            self.dimension = len(self.coordinates)

        except ValueError:
            raise ValueError('The coordinates must be nonempty')

        except TypeError:
            raise TypeError('The coordinates must be an iterable')


    def __str__(self):
        return 'Vector: {}'.format(self.coordinates)


    def __eq__(self, v):
        return self.coordinates == v.coordinates
        
    def plus(self, v):
        new_coordinates = [x+y for x,y in zip(self.coordinates,v.coordinates)]
        return Vector(new_coordinates)
        
    def minus(self, v):
        new_coordinates = [x-y for x,y in zip(self.coordinates,v.coordinates)]
        return Vector(new_coordinates)
    
    def times_scalar(self, c):
        new_coordinates = [Decimal(c)*x for x in self.coordinates]
        return Vector(new_coordinates)
        
    def magnitude(self):
        return (sum([v**Decimal('2') for v in self.coordinates])).sqrt()
        
    def normalize(self):
        try:
            magnitude = self.magnitude()
            return self.times_scalar(Decimal('1.0')/magnitude)
            
        except ZeroDivisionError:
            raise Exception('Cannot normalize the zero vector')
            
    def dot_product(self,v):
        dot = sum([x*y for x,y in zip(self.coordinates,v.coordinates)])
        return dot
        
    def theta(self,v,in_degrees=False):
        try:
            mag_self, mag_v = self.magnitude(), v.magnitude()
            angle = acos(self.dot_product(v) / (mag_self * mag_v))
            
            if in_degrees:
                degrees_per_radian = 180./pi
                return angle * degrees_per_radian
            else:
                return angle
            
        except ZeroDivisionError:
            raise Exception('Cannot divide by zero')
            
    def is_parallel(self,v):
        check_coordinates = {x/y if y != 0 else 0 for x,y in zip(self.coordinates,v.coordinates)}
        return len(check_coordinates) == 1
        
    def is_ortho(self,v, tolerance=1e-10):
        return abs(self.dot_product(v)) < tolerance    

    def projection(self,b):
        try:
            unit_vector_b = b.normalize()
            x = self.dot_product(unit_vector_b)
            return unit_vector_b.times_scalar(x)
        except Exception as e:
            if str(e) == 'Cannot normalize the zero vector':
                raise Exception('No unique parallel component to zero vector')
            else:
                raise e
        
        
    def ortho_component(self,b):
        v_proj = self.projection(b)
        return self.minus(v_proj)
        
    def decompose(self,b):
        proj = self.projection(b)
        comp = self.ortho_component(b)
        return (proj,comp)
        
        
        
        