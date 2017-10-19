__precompile__() 

module herramientas

using SymPy

v=Sym("v");
w=Sym("w");

export Newton
export Newton2
export riemann
export trapecio
export simpson
export lagrange
export numDer
export euler1D
export eulerND
export impEuler
export midP
export RK
export RK4

##Se omitieron todos los comentarios, pues todas las funciones ya estaban definidas en otros notebooks

"Método de Newton para encontrar raíces."""
#Argumentos (función, derivada de la función, suposición incial)

function Newton(f, df, y)
       x=y
       for i in 1:200
       x = x-(f(x)/df(x))
       end
       return x
       end
       
"""Método de Newton actualizado para que sirva con las funciones F y dF del método implícito de Euler.""" #Argumentos(función en la que se aplica el método, su derivada, la suposición de raíz, t es el tiempo para la función de dos variables)

function Newton2(F, dF, y, t, y0, h)
    x=y
       for i in 1:2000
       x = x-(F(x, y0, h, t)/dF(x, y0, h, t))
       end
       return x;
end;

"""Integración con Riemann."""
#Argumentos (función, x inicial, x final, tamaño del subintervalo)

function riemann(f::Function, a, b, delta)
    integral = 0;
   
    n = round((b-a)/delta);
   
    an=a;
    bn=an+delta;
    
    for i in 1:n
        integral+=delta*f((an+bn)/2);
        an = bn;
        bn = bn+delta;
    end
    
    return integral;
    
end

"""Integración con método del trapecio."""
#Argumentos (función, x inicial, x final, tamaño del subintervalo)

function trapecio(f::Function, a, b, delta)
    integral = 0;
    
    n = (b-a)/delta;
   
    an=a;
    bn=an+delta;
    
    
    for i in 1:n
        integral+=delta*((f(an)+f(bn))/2);
        an = bn;
        bn = bn+delta;
    end
    
    return integral;
    
end

"""Integración con Simpson."""
#Argumentos (función, x inicial, x final, tamaño del subintervalo)

function simpson(f::Function, a, b, delta)
    integral = 0;
    
    n = (b-a)/delta;
    
    an=a;
    bn=an+delta;
    
   
    for i in 1:n
        integral+=(delta/6)*(f(an)+4*f((an+bn)/2)+f(bn));
        an = bn;
        bn = bn+delta;
    end
    
    return integral;
    
end

"""Interpolación de Lagrange, regresa el polinomio."""
#Argumentos(arreglo de puntos en x, arreglo de puntos en y, punto en donde se evalúa el polinomio, puede ser variable simbólica)

function lagrange(a, b, eval)
  
    xArray = collect(a);
    yArray = collect(b);
    
  
    x = Sym("x");
    
   
    L=0*x;
    
    
    for j in 1:length(xArray)
        n=1
       
        h=1.0*x/x;
       
        while n<=length(xArray)
           
           if n==j
           n = n+1;
           else
         
            h= h*(x-xArray[n])/(xArray[j]-xArray[n]);
            n = n+1;
           end
        end
        
        L = L + yArray[j]*h;
    end
    
  
    L = simplify(L);
    
   
    pol=lambdify(L,[x]);
    
    
    return pol(eval);
    
end

"""Derivada numérica, regresa el valor en el punto dado.""" 
#Argumentos(función, punto en donde se evalúa, tamaño de h)

function numDer(f::Function, x, h)
  
    temp = (f(x+h)-f(x))/h;
    return temp;
end


"""Método de Euler en 1D, regresa dos arrays el primero de tiempo y el segundo de posición."""
#Argumentos(función, condición inicial, primer tiempo, último tiempo y el tamaño de cada subintervalo)

function euler1D(g, a, ti, tf, h)
    
    timeArray = linspace(ti, tf, round((tf-ti)/h));
    xArray = zeros(length(timeArray));
   
    xArray[1] = a;

    for i in 1:length(timeArray)-1
        xArray[i+1] = xArray[i] +h*g(timeArray[i], xArray[i]);
    end
   
    return timeArray, xArray; 
end;

"""Método de Euler para varias dimensiones, regresa una lista de posición."""
#Argumentos(función, lista de tiempo, condición inicial)

function eulerND(f,list,x0)
     x = x0
     h = list[2]-list[1]
     listx=[]
     push!(listx,x)
     for i in 2:length(list)
        t = i*h
        x = x + f(x,t)*h
        push!(listx,x) 
     end
     return listx
end


"""Método de Euler implícito, regresa dos arrays el primero de tiempo y el segundo de posición."""
#Argumentos(función, derivada de la función, tiempo inicial, tiempo final, condición inicial, tamaño de subintervalo) Recibe la función F=y[n+1]-y[n]-h*f(y[n],t[n])

function impEuler(F, dF, t0, tn, y0, h)
    
    n = trunc(Int,(tn - t0)/h);
   
    timeArray = zeros(n+1);
    yArray = zeros(n+1);
    
    timeArray[1] = t0
    yArray[1] = y0
   
    for i in 1:n
        timeArray[i+1]=timeArray[i] + h;
        yArray[i+1] = rootNewtonMethod2(F,dF,1,timeArray[i+1], yArray[i], h);
    end;
    
    return timeArray, yArray;
end;

"""Regla del punto medio, regresa dos arrays el primero de tiempo y el segundo de posición. """
#Argumentos(función, condición incial, tiempo inicial, tiempo final, tamaño de subintervalo)

function midP(g, a, ti, tf, h)
   
    timeArray = linspace(ti, tf, round((tf-ti)/(h)));
    xArray = zeros(length(timeArray));
   
    xArray[1] = a;
    
    for i in 1:length(timeArray)-1
        xArray[i+1] = xArray[i] +h*g(xArray[i]+(h/2)*g(xArray[i],timeArray[i]), timeArray[i]+(h/2));
    end
   
    return timeArray, xArray; 
end;

"""Runge-Kutta de orden 4 para 1D, regresa dos arrays el primero de tiempo y el segundo de posición."""
#Argumentos(función, condición incial, tiempo inicial, tiempo final, tamaño de subintervalo)

function RK(g, a, ti, tf, h)
    
    timeArray = linspace(ti, tf, round((tf-ti)/(h)));
    xArray = zeros(length(timeArray));
    
    xArray[1] = a;
    
    for i in 1:length(timeArray)-1
        k1 = g(xArray[i], timeArray[i]);
        k2 = g(xArray[i]+(h/2)*k1, timeArray[i]+(h/2));
        k3 = g(xArray[i]+(h/2)*k2, timeArray[i]+(h/2));
        k4 = g(xArray[i]+h*k3, timeArray[i+1]);
        xArray[i+1] = xArray[i]+(h/6)*(k1+2*k2+2*k3+k4);
    end
    
    return timeArray, xArray; 
end;

"""Runge-Kutta de orden 4 para cualquier dimensión, regresa una lista de posición.""" 
#Argumentos(función, lista de tiempo, condición inicial)

function RK4(f,list,x0)
     x = x0
     h = list[2]-list[1]
     listx=[]
     push!(listx,x)
     for i in 2:length(list)
        t = i*h
        k1 = f(x,t);
        k2 = f(x+(h/2)*k1,t+(h/2));
        k3 = f(x+(h/2)*k2, t+(h/2));
        k4 = f(x+h*k3, t+h);
        x = x+(h/6)*(k1+2*k2+2*k3+k4);
        push!(listx,x) 
     end
     return listx
end



end