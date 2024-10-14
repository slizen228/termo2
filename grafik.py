import numpy as np
import plotly.graph_objs as go

N = int(input("N"))+1
M = int(input("M"))+1


# Допустим, у вас уже есть массив значений функции u (размер 10x10)
u1 = np.loadtxt('value.txt')
u2 = np.loadtxt('value_function.txt')


# Создаем массивы для координат
i = np.arange(N)
j = np.arange(M)

# Преобразуем координаты в реальные значения
x, y = np.meshgrid(i * (1/(N-1)), j * (1/(M-1)))

# Создаем 3D поверхность для отображения
surface1 = go.Surface(z=u1, x=x, y=y, colorscale='blues', name='найденное решение')

surface2 = go.Surface(z=u2, x=x, y=y, colorscale='Cividis', name='точное решение')

# Настраиваем макет с явными осями координат
layout = go.Layout(
    title='3D график с двумя функциями',
    scene=dict(
        xaxis=dict(title='X (i *(1/(N-1)))', showgrid=True, zeroline=True, showline=True, ticks='outside', showticklabels=True),
        yaxis=dict(title='Y (j * (1/(M-1))', showgrid=True, zeroline=True, showline=True, ticks='outside', showticklabels=True),
        zaxis=dict(title='U (значения функций)', showgrid=True, zeroline=True, showline=True, ticks='outside', showticklabels=True)
    )
)

# Создаем фигуру с двумя поверхностями
fig = go.Figure(data=[surface1,surface2], layout=layout)

# Визуализируем график
fig.show()
