import math
import time
import ovito
from ovito.io import import_file
from ovito.modifiers import ConstructSurfaceModifier
from ovito.modifiers import DislocationAnalysisModifier, SelectExpressionModifier, DeleteSelectedParticlesModifier

print("Disloc script is started")

# Extract dislocation lines from the crystal:
st="/mnt/disk1/GB/Erate/results/Al_GB_5_19_421_erate.dump"

node = import_file(st, multiple_frames = True)
#node = import_file(st, multiple_frames = False)

g = open('gistogramma_dislocations_10_38_421.txt', 'a')

#g.write("step time koncentraziya[0] koncentraziya[1] koncentraziya[2] Shockley Perfect Stair Frank Hirth rhoD RHOD n_seg\n")
#g.write("step time koncentraziya[0] koncentraziya[1] koncentraziya[2] koncentraziya[3] koncentraziya[4] koncentraziya[5] koncentraziya[6] koncentraziya[7] koncentraziya[8] koncentraziya[9] n_seg\n")
g.write("step time layer_x RHOD\n")

mod = ConstructSurfaceModifier(radius = 2.4)
modDXA = DislocationAnalysisModifier()
node.modifiers.append(mod)
node.modifiers.append(modDXA)

modifier =  DislocationAnalysisModifier()
modifier.input_crystal_structure = DislocationAnalysisModifier.Lattice.FCC
node.modifiers.append(modifier)

def function_peresecheniya(x1, y1, z1, x2, y2, z2, x_pl1, y_pl1, z_pl1, x_pl2, y_pl2, z_pl2, x_pl3, y_pl3, z_pl3):
   #векторы лежащие в плоскости
   #компоненты первого вектора
   ax = x_pl2 - x_pl1
   ay = y_pl2 - y_pl1
   az = z_pl2 - z_pl1
   #компоненты второго вектора
   bx = x_pl3 - x_pl1
   by = y_pl3 - y_pl1
   bz = z_pl3 - z_pl1
   #компоненты вектора нормали к плоскости
   i = ay*bz-az*by 
   j = -(ax*bz-az*bx)
   k = ax*by-ay*bx 
   #уравнение плоскости  Ax+By+Cz+D=0
   A = i
   B = j
   C = k
   D = i * x_pl1 + j * y_pl1 + k * z_pl1
   #параметрическое уравнение прямой   x=x0+m*t    y=y0+n*t    z=z0+p*t
   m = x2 - x1
   n = y2 - y1
   p = z2 - z1
   t = (A * x1 + B * y1 + C * z1 + D)/(A * m + B * n + C * p)   
   #точки пересечения прямой с плоскостью 
   x=x1+m*t
   y=y1+n*t
   z=z1+p*t  
   dlina_chasti_vectora = math.sqrt((x - x1)**2 + (y - y1)**2 + (z - z1)**2)
   return dlina_chasti_vectora

f_start = 69

for f in range(f_start, ovito.dataset.anim.last_frame+1):
    ovito.dataset.anim.current_frame = f
    node.compute()
    print("frame=", f)
    time = f * 5.0

    # Print list of dislocation lines:
    network = node.output.dislocations

    print("Found %i dislocation segments " % len(network.segments))

    # DXA
    #solid_volume = node.output.attributes['ConstructSurfaceMesh.solid_volume']
    total_line_length = node.output.attributes['DislocationAnalysis.total_line_length']/10.0
    number_sigments = len(node.output.dislocations.segments)
    counts_HCP = node.output.attributes['DislocationAnalysis.counts.HCP']

    vol=node.output.attributes['DislocationAnalysis.cell_volume']
    Shockley=node.output.attributes['DislocationAnalysis.length.1/6<112>'] / vol
    Perfect=node.output.attributes['DislocationAnalysis.length.1/2<110>'] / vol
    Stair=node.output.attributes['DislocationAnalysis.length.1/6<110>'] / vol
    Frank=node.output.attributes['DislocationAnalysis.length.1/3<111>'] / vol
    Hirth=node.output.attributes['DislocationAnalysis.length.1/3<100>'] / vol

    #RHOD = total_line_length*10.0 / solid_volume * 1.0e+8
    rhoD = 0


    #размер системы по X
    length_system_x = 1042
    length_system_z = 359
    #на сколько слоев разбиваем систему по x
    number_layer = 50

    #координаты трех точек образующих плоскость
    y_plane1 = 20
    z_plane1 = 30
    y_plane2 = 120
    z_plane2 = 300
    y_plane3 = 200
    z_plane3 = 170

    #задание размеров одного слоя
    x_length_sloya = length_system_x/number_layer
    y_length_sloya = 230
    z_length_sloya = 359
    volume = x_length_sloya * y_length_sloya * z_length_sloya

    #вычисление дельты и угла для порпавки на смещение
    velocity_deformation = 0.000327
    delta_x =  length_system_z*velocity_deformation*1e12 * time * 1e-12
    Alpha =  math.acos(delta_x / math.sqrt(delta_x**2 + z_length_sloya**2))


    koncentraciya = [0] * number_layer
    RHOD = [0] * number_layer
    length_segment = 0
    #пробегаем по сегментам в системе
    for segment in network.segments:
        x_beg = 0
        y_beg = 0
        z_beg = 0
        x_end = 0
        y_end = 0
        z_end = 0
        #проходим по каждой точке в одном сегменте
        for point in segment.points:
            #задаем координаты для второй точки при второй и следующих итерациях
            x_beg = point[0]
            y_beg = point[1]
            z_beg = point[2]
            #print("x_beg=",x_beg)
            #print("y_beg=",y_beg)
            #print("z_beg=",z_beg)
            #print("x_end=",x_end)
            #print("y_end=",y_end)
            #print("z_end=",z_end)
            #проверяем чтобы первая точка была ненулевой и если выполняется условие можно посчитать длину
            if x_end != 0 or y_end != 0 or z_end != 0:
               #пробегаем по слоям системы
               #print("зашли в условие")
               counter = 0
               for i in range(-1, number_layer, 1):
                  #вычисляем смещение изза сдвиговой деформации
                  x_fix1 = z_end / math.tan(Alpha)
                  x_fix2 = z_beg / math.tan(Alpha)
                  dlina_megdu_tochkami = 0
                  if i < 0:
                     counter = 0
                  elif i >= number_layer:
                     counter =  number_layer - 1
                  else:
                     counter = i
                  #print(i)
                  #ставим условие если при итерировании слоев обе точки попадают в один слой то можно посчитать длину и записать в массив
                  if x_end > (i*x_length_sloya) + x_fix1 and x_end < (i+1)*x_length_sloya + x_fix1 and x_beg> i*x_length_sloya + x_fix2 and x_beg < (i+1)*x_length_sloya + x_fix2:
                     #print("зашли в условие 1")
                     #print("координата точки x_beg=",x_beg)
                     #print("координата точки x_end=",x_end)

                     dlina_megdu_tochkami = math.sqrt((x_end - x_beg)**2 + (y_end - y_beg)**2 + (z_end - z_beg)**2)
                     #print("длина между точками = ",dlina_megdu_tochkami)
                     #print("заносим в koncentraciya["+str(i)+"] =")

                     koncentraciya[counter] += dlina_megdu_tochkami
                     #print("заносим в koncentraciya["+str(i)+"] =" , koncentraciya[i])
                     #print("заносим точки в слой["+str(i)+"]")
                     #time.sleep(5)
                     #print("поправка по z координате равна=",x_fix1 )
                     #print("поправка по z координате равна=",x_fix2 )
                  #иначе находим в каких слоях находятся точки
                  elif x_end > i*x_length_sloya + x_fix1 and x_end < (i+1)*x_length_sloya + x_fix1:  #находим первую точку
                     #запускаем цикл для нахождения слоя в котором присутствует вторая точка
                     for j in range (-1, number_layer, 1):
                        if x_beg> (j*x_length_sloya) + x_fix2 and x_beg < ((j+1)*x_length_sloya) + x_fix2:
                           if j < 0:
                              iterator = 0
                           elif j >= number_layer:
                              iterator = number_layer - 1
                           else:
                              iterator = j
                           n = iterator - i  #kolichestvo ploskostey megdu tochkami
                           # x координаты точек которые образуют текущую плоскоть
                           x_plane1 = 0
                           x_plane2 = 0
                           x_plane3 = 0
                           #print("зашли в условие 1")
                           #print("координата точки x_beg=",x_beg)
                           #print("координата точки x_end=",x_end)
                           #проверяем с какой стороны находится вторая точка относительно первой
                           if n < 0: # koordinata ploskosti sleva
                              #запускаем цикл и итерируем слои которые находятся между точками
                              dlina_ostatochnoy_chasty = 0
                              for l in range(counter, iterator-1, -1):
                                 # koordinata ploskosti sleva
                                 x_plane1 = l*x_length_sloya + z_plane1 / math.tan(Alpha)
                                 x_plane2 = l*x_length_sloya + z_plane2 / math.tan(Alpha)
                                 x_plane3 = l*x_length_sloya + z_plane3 / math.tan(Alpha)
                                 #вычисляем длину вектора вплоть до пересечения ее с плоскостью
                                 dlina_vectora = function_peresecheniya(x_end, y_end, z_end, x_beg, y_beg, z_beg, x_plane1, y_plane1, z_plane1, x_plane2, y_plane2, z_plane2, x_plane3, y_plane3, z_plane3)
                                 #записываем эту часть длины вектора в определенный слой
                                 koncentraciya[l] = koncentraciya[l] + (dlina_vectora - dlina_ostatochnoy_chasty)
                                 #print("koordinata ploskosti sleva ")
                                 #print("длина между точками = ",dlina_vectora)
                                 #print("заносим в koncentraciya["+str(l)+"] =")
                                 #print("заносим в koncentraciya["+str(l)+"] =" , koncentraciya[l])
                                 #записываем остаточную часть вектора чтобы убрать кусочек длины уже записанный в предыдущий слой
                                 dlina_ostatochnoy_chasty = dlina_vectora
                                 #print("заносим точки в слой["+str(l)+"]")
                           elif n > 0: # koordinata ploskosti sprava
                              #запускаем цикл и итерируем слои которые находятся между точками
                              dlina_ostatochnoy_chasty = 0
                              for l in range(counter, iterator+1, 1):
                                 # koordinata ploskosti sprava
                                 x_plane1 = (l+1)*x_length_sloya + z_plane1 / math.tan(Alpha)
                                 x_plane2 = (l+1)*x_length_sloya + z_plane2 / math.tan(Alpha)
                                 x_plane3 = (l+1)*x_length_sloya + z_plane3 / math.tan(Alpha)
                                 #вычисляем длину вектора вплоть до пересечения ее с плоскостью
                                 dlina_vectora = function_peresecheniya(x_end, y_end, z_end, x_beg, y_beg, z_beg, x_plane1, y_plane1, z_plane1, x_plane2, y_plane2, z_plane2, x_plane3, y_plane3, z_plane3)
                                 #записываем эту часть длины вектора в определенный слой
                                 koncentraciya[l] = koncentraciya[l] + (dlina_vectora - dlina_ostatochnoy_chasty)
                                 #print("koordinata ploskosti sprava ")
                                 #print("длина между точками = ",dlina_vectora)
                                 #print("заносим в koncentraciya["+str(l)+"] =")
                                 #print("заносим в koncentraciya["+str(l)+"] =" , koncentraciya[l])
                                 #записываем остаточную часть вектора чтобы убрать кусочек длины уже записанный в предыдущий слой
                                 dlina_ostatochnoy_chasty = dlina_vectora
                                 #print("заносим точки в слой["+str(l)+"]")

            #задаем координаты для первой точки
            x_end = point[0]
            y_end = point[1]
            z_end = point[2]
    print("\n frame= ", f)
    print("delta_x= ", delta_x)
    print("Alpha= ", Alpha)
    print("\n")
    for k in range(0, number_layer, 1):
       RHOD[k] = koncentraciya[k]/volume
       print("RHOD"+str(k)+"=", RHOD[k])
       g.write("%i %f %i %f  \n" % (f, time, k, RHOD[k]))

    #print("step=%i x_a=%f y_a=%f z_a=%f n=%f" % (f, x_aver, y_aver, z_aver, n))
    #g.write("%i %f %i %i %i %f %f %f %f %f %f %f %i\n" % (f, time, koncentraziya[0], koncentraziya[1], koncentraziya[2], Shockley, Perfect, Stair, Frank, Hirth, rhoD, RHOD, n))

print ("Done")
g.close()
