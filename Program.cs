using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Threading.Tasks;
using Tao.OpenGl;
using Tao.FreeGlut;
using System.Threading;

using System.Windows.Input;
using System.Windows;
using System.Collections;


//1 - yellow, 2- red, 3- blue, 4- green, else - black
namespace Graph_Console
{
    class Program
    {
        static Graph A = new Graph();


        const int NMax = 20;
        const double R_bigCircle = 0.7;
        const double r_smallCircle = 0.07;
        const double h_num = 0.018;
        const int INF = 1000000;

        public class edge
        {
            public int from;
            public int to;
            public int weight;
           // public int color;

        }

        public class vertex
        {
            public List<edge> rebra = new List<edge>();
            public int number;
            public int color;
            public edge MinRebro()
            {
                edge min = new edge();
                min = rebra[0];
                for (int j = 0; j < rebra.Count; j++)
                    if (min.weight > rebra[j].weight)
                        min = rebra[j];

                return min;
            }



        }

        public class Graph
        {
            public List<edge> seq_e = new List<edge>();
            public List<int> seq = new List<int>();
           // public int currentStep1, currentStep2;
           // public DateTime lastStepDate;
            //public bool ShowAnimation;


            public int[,] matrix;
            public int[,] color_Edge;//2
            public int[] color_Vertex;//1

            public void Initialize()
            {

                ReadFromFileLineByLineIntoArray(ref matrix);
                color_Edge = new int[matrix.GetLength(0), matrix.GetLength(0)];
                color_Vertex = new int[matrix.GetLength(0)];
                for (int i = 0; i < matrix.GetLength(0); i++) color_Vertex[i] = 1;
                for (int j = 0; j < matrix.GetLength(0); j++)
                    for (int k = 0; k < matrix.GetLength(1); k++)
                    {
                        if (matrix[j, k] != 0) color_Edge[j, k] = 2;
                        else matrix[j, k] = INF;
                    }



                Console.Write(" \nVertex\n");
                foreach (var element in color_Vertex) Console.Write(element + " ");
                Console.Write(" \nEdge\n");

                for (int k = 0; k < color_Edge.GetLength(0); ++k)
                {
                    for (int j = 0; j < color_Edge.GetLength(1); ++j)
                    { Console.Write(color_Edge[k, j] + " "); }
                    Console.WriteLine();
                }

                Console.WriteLine();
            }
            public void First()
            {
                for (int i = 0; i < matrix.GetLength(0); i++) color_Vertex[i] = 1;
                for (int j = 0; j < matrix.GetLength(0); j++)
                    for (int k = 0; k < matrix.GetLength(1); k++)
                        if (matrix[j, k] != INF) color_Edge[j, k] = 2;
            }

 
 
            public void Draw()
            {

                Gl.glClear(Gl.GL_COLOR_BUFFER_BIT | Gl.GL_DEPTH_BUFFER_BIT);

                int count = matrix.GetLength(0);
                double psi = (2 * Math.PI) / count;
                bool or = false;
                for (int i = 0; i < matrix.GetLength(0); ++i)
                    for (int j = i; j < matrix.GetLength(1); ++j)//symmetric
                    { if (matrix[i, j] != matrix[j, i]) { or = true; break; } }

                for (int i = 0; i < matrix.GetLength(0); ++i)
                    for (int j = i; j < matrix.GetLength(1); ++j)//symmetric
                        if (matrix[i, j] != INF)
                        {
                            double x1 = R_bigCircle * Math.Cos(psi * i);
                            double y1 = R_bigCircle * Math.Sin(psi * i);
                            double x2 = R_bigCircle * Math.Cos(psi * (j));
                            double y2 = R_bigCircle * Math.Sin(psi * (j));
                            double a = (x2 - x1);
                            double b = (y2 - y1);
                            double d = Math.Sqrt((y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1));
                            double ex = a / d;
                            double ey = b / d;
                            Picture.Line(x1, y1, x2, y2, color_Edge[i, j]);
                            if (or)
                            {
                                double x3 = x2- ex * r_smallCircle;
                                double y3 = y2 - ey * r_smallCircle;
                                double d3 = 0.02;
                                double nx1 = y2 - y1;
                                double ny1 = x1 - x2;
                                double vecy = y1 - y2;
                                double vecx = x1 - x2;
                                double nd = Math.Sqrt((ny1 * ny1) + (nx1 * nx1));
                                double vecd = Math.Sqrt((vecy * vecy) + (vecx * vecx));

                                Picture.Line(x3, y3, x3 + (vecx / vecd) * d3 + (nx1 / nd) * d3, y3 + (vecy / vecd) * d3 + (ny1 / nd) * d3, color_Edge[i, j]);
                                Picture.Line(x3, y3, x3 + (vecx / vecd) * d3 - (nx1 / nd) * d3, y3 + (vecy / vecd) * d3 - (ny1 / nd) * d3, color_Edge[i, j]);
                              
                            }
                            Picture.DrawNumbersInverse(matrix[i, j], ex * (d / 2.0) + x1 + 0.02, ey * (d / 2.0) + y1 + 0.02,0);
                        }
                for (int i = 0; i < count; i++)
                {
                    Picture.Circle(R_bigCircle * Math.Cos(psi * i), R_bigCircle * Math.Sin(psi * i), color_Vertex[i]);
                    Picture.DrawNumbersInverse(i, R_bigCircle * Math.Cos(psi * i), R_bigCircle * Math.Sin(psi * i),0);
                }
                 // Gl.glFlush();
                Glut.glutSwapBuffers();
            }
            public void Draw_flow( int[,] capacity, int [,] flow)
            {

                Gl.glClear(Gl.GL_COLOR_BUFFER_BIT | Gl.GL_DEPTH_BUFFER_BIT);

                int count = matrix.GetLength(0);
                double psi = (2 * Math.PI) / count;

                bool or = false;
                for (int i = 0; i < matrix.GetLength(0); ++i)
                    for (int j = i; j < matrix.GetLength(1); ++j)
                    { if (matrix[i, j] != matrix[j, i]) { or = true; break; } }

                for (int i = 0; i < matrix.GetLength(0); ++i)
                    for (int j = 0; j < matrix.GetLength(1); ++j)
                        if (matrix[i, j] != INF)
                        {
                            double x1 = R_bigCircle * Math.Cos(psi * i);
                            double y1 = R_bigCircle * Math.Sin(psi * i);
                            double x2 = R_bigCircle * Math.Cos(psi * (j));
                            double y2 = R_bigCircle * Math.Sin(psi * (j));
                            double a = (x2 - x1);
                            double b = (y2 - y1);
                            double d = Math.Sqrt((y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1));
                            double ex = a / d;
                            double ey = b / d;
                            Picture.Line(x1, y1, x2, y2, color_Edge[i, j]);
                            Picture.DrawNumbersInverse(capacity[i,j], ex * (d / 2.0) + x1 + 0.02, ey * (d / 2.0) + y1 + 0.02, 0);
                            Picture.DrawNumbersInverse(flow[i, j], ex * (d / 2.0) + x1 - 0.02, ey * (d / 2.0) + y1 - 0.02, 3);

                            if (or)
                            {
                                double x3 = x2 - ex * r_smallCircle;
                                double y3 = y2 - ey * r_smallCircle;
                                double d3 = 0.02;
                                double nx1 = y2 - y1;
                                double ny1 = x1 - x2;
                                double vecy = y1 - y2;
                                double vecx = x1 - x2;
                                double nd = Math.Sqrt((ny1 * ny1) + (nx1 * nx1));
                                double vecd = Math.Sqrt((vecy * vecy) + (vecx * vecx));

                                Picture.Line(x3, y3, x3 + (vecx / vecd) * d3 + (nx1 / nd) * d3, y3 + (vecy / vecd) * d3 + (ny1 / nd) * d3, color_Edge[i, j]);
                                Picture.Line(x3, y3, x3 + (vecx / vecd) * d3 - (nx1 / nd) * d3, y3 + (vecy / vecd) * d3 - (ny1 / nd) * d3, color_Edge[i, j]);

                            }
                        }
                for (int i = 0; i < count; i++)
                {
                    Picture.Circle(R_bigCircle * Math.Cos(psi * i), R_bigCircle * Math.Sin(psi * i), color_Vertex[i]);
                    Picture.DrawNumbersInverse(i, R_bigCircle * Math.Cos(psi * i), R_bigCircle * Math.Sin(psi * i), 0);
                }
                // Gl.glFlush();
                Glut.glutSwapBuffers();
            }

            public void DFS(int number)  // после дфс обязательно консоль ридлайн иначе не видно последней вершины( рек)
            {

                seq.Add(number);

               Console.ReadLine();
                color_Vertex[number] = 3;//blue
                Draw();

                for (int j = 0; j < matrix.GetLength(0); j++)
                    if (matrix[number, j] != INF)
                        if (color_Vertex[j] == 1) DFS(j);

            }
            public bool Cicle(int number)
            {

                //seq.Add(number);

                color_Vertex[number] = 3;//blue


                for (int j = 0; j < matrix.GetLength(0); j++)
                {
                    if (matrix[number, j] != INF)
                    {
                        int to = j;
                        if (color_Vertex[to] == 1)
                        {
                            seq.Add(number);
                            if (Cicle(to)) return true;
                        }



                        else
                        {
                            if (((color_Vertex[to] == 3) && (to != seq[seq.Count - 1])))
                                return true;
                        }
                    }
                }
                //   color_Vertex[number] = 2;
                return false;
            }
            public void BFS(int number_v)
            {


                color_Vertex[number_v] = 3;
                Draw();
                Console.ReadLine();
                int t = 0, h = 0;
                int[] queue = new int[matrix.GetLength(0)];
                queue[t] = number_v; t++;
                while (h < t)
                {
                    int v = queue[h]; h++;
                    for (int j = 0; j < matrix.GetLength(0); ++j)
                    {
                        if (matrix[v, j] != INF)
                        {
                            if (color_Vertex[j] == 1)
                            {
                                color_Vertex[j] = 3;
                                queue[t] = j; t++;
                                Draw();
                                Console.ReadLine();

                            }
                        }
                    }
                }
                for (int i = 0; i < t; i++)
                {
                    seq.Add(queue[i]);
                }

            }
            public void Kruskala()
            {
                List<edge> wEdge;
                List<edge> temp = new List<edge>();
                CreateListEdge(out wEdge);
                SortByWeightListEdge(wEdge);
                int index = 0;
                 seq.Clear();
                for (int i = 0; i < wEdge.Count; i++)
                {

                    Graph gr = new Graph();
                    temp.Add(wEdge[i]);
                    gr.FromListOfEdgesToGraph(temp);
                    if (!gr.Cicle(0))
                    {
                        color_Edge[wEdge[i].from, wEdge[i].to] = 4;
                        color_Edge[wEdge[i].to, wEdge[i].from] = 4;

                    }
                    else { temp.RemoveAt(index); index--; }
                    index++;
                     Draw();
                     Console.ReadLine();
                    
                }
                seq_e.AddRange(temp);
            }
            public void Prima(int num_v)
            {

                
                int[] min_e = new int[matrix.GetLength(0)];
                int[] sel_e = new int[matrix.GetLength(0)];
                bool[] used = new bool[matrix.GetLength(0)];

                for (int i = 0; i < matrix.GetLength(0); i++)
                {
                    min_e[i] = INF;
                    sel_e[i] = -1;
                    used[i] = false;
                }

                min_e[num_v] = 0;
                for (int i = 0; i < matrix.GetLength(0); ++i)
                {
                    int v = -1;
                    for (int j = 0; j < matrix.GetLength(0); ++j)
                        if (!used[j] && (v == -1 || min_e[j] < min_e[v])) v = j;
                    if (min_e[v] == INF) { Console.WriteLine("No MST!"); return; }

                    used[v] = true;
                    color_Vertex[v] = 3;
                    seq.Add(v);

                    if (sel_e[v] != -1)
                    {
                        Console.WriteLine("sel_e[{0}] = {1} ", v, sel_e[v]);
                        color_Edge[sel_e[v], v] = 4;
                        color_Edge[v, sel_e[v]] = 4;
                        edge w = new edge();
                        w.from = sel_e[v];
                        w.to = v;
                       // w.color = 4;
                        w.weight = matrix[sel_e[v], v];
                        seq_e.Add(w);
                    }
                     Draw();
                     Console.ReadLine();

                    for (int to = 0; to < matrix.GetLength(0); ++to)
                        if (matrix[v, to] < min_e[to])
                        {
                            min_e[to] = matrix[v, to];
                            sel_e[to] = v;
                        }
                }
            }
            public void Deikstra(int num_v)
                {
                    int[] d = new int[matrix.GetLength(0)];
                    int[] p = new int[matrix.GetLength(0)];
                    int[] sel_e = new int[matrix.GetLength(0)];
                    bool[] used = new bool[matrix.GetLength(0)];

                    for (int i = 0; i < matrix.GetLength(0); i++)
                    {
                        d[i] = INF;
                        sel_e[i] = -1;
                        p[i]=-1;
                        used[i] = false;
                    }


                    d[num_v] = 0;

                    for (int i = 0; i < matrix.GetLength(0); ++i)
                    {
                        int v = -1;
                        for (int j = 0; j < matrix.GetLength(0); ++j)
                            if (!used[j] && (v == -1 || d[j] < d[v]))
                                v = j;
                        if (d[v] == INF) { Console.WriteLine("No MST!"); return; }

                        used[v] = true;
                        color_Vertex[v] = 3;


                        if (p[v] != -1)
                        {
                            Console.WriteLine("p[{0}] = {1} ", v, p[v]);

                        }
                        Draw();
                        Console.ReadLine();

                        for (int k = 0; k < matrix.GetLength(0); ++k)
                        {
                            if ((matrix[v, k] != INF) && (!used[k]))
                            {
                                int to = k;
                                int len = matrix[v, k];
                                if (d[v] + len < d[to])
                                {
                                    d[to] = d[v] + len;
                                    p[to] = v;
                                }
                            }
                        }
                    }
                        double psi = (2 * Math.PI) / matrix.GetLength(0);
                        for (int i = 0; i < matrix.GetLength(0); i++)
                        {
                            
                            Picture.DrawNumbersInverse(d[i], R_bigCircle * Math.Cos(psi * i)+r_smallCircle, R_bigCircle * Math.Sin(psi * i)+r_smallCircle,6);
                        }
                        Gl.glFlush();
                        Console.ReadLine();

                    
            }
            public void FloydWarshall()//or
            {
                int[,] W = new int[matrix.GetLength(0), matrix.GetLength(0)];
                for (int i = 0; i < matrix.GetLength(0); i++)
                    for (int j = 0; j < matrix.GetLength(1); j++)
                    {
                       
                        W[i, j] = matrix[i, j];
                        
                    }
                    

                for(int k=0; k< matrix.GetLength(0); k++)
                    for( int i=0; i< matrix.GetLength(0); i++)
                        for (int j = 0; j < matrix.GetLength(0); j++)  
                            
                            W[i,j] = min(W[i,j], W[i,k] + W[k,j]);
             Console.WriteLine("FloydWarshall algorithm\n\n");
             for (int k = 0; k < matrix.GetLength(0); ++k)
             {
                 for (int j = 0; j < matrix.GetLength(1); ++j)
                 { Console.Write(W[k,j] + " "); }
                 Console.WriteLine();
             }

             Console.WriteLine();
            
            }
            public void BellmanFord(int s)//or
            {
                int[] d = new int[matrix.GetLength(0)];
                for (int i = 0; i < d.Length; i++)
                    d[i] = INF;

                d[s] = 0;
                int l=0;
                for (; l < matrix.GetLength(0); l++)
                {
                    for (int j = 0; j<matrix.GetLength(0); j++)
                        for( int k=0; k<matrix.GetLength(1); k++)
                            if (matrix[j, k] != INF)
                            {
                                if (d[k] > d[j] + matrix[j, k]) d[k] = d[j] + matrix[j, k];
                            }
                }
                int [] check= new int [matrix.GetLength(0)];
                for (int i = 0; i < d.Length; i++)
                {
                    check[i] = d[i];
                }
                for (; l < matrix.GetLength(0)+1; l++)
                {
                    for (int j = 0; j < matrix.GetLength(0); j++)
                        for (int k = 0; k < matrix.GetLength(1); k++)
                            if (matrix[j, k] != INF)
                            {
                                if (check[k] > check[j] + matrix[j, k]) check[k] = check[j] + matrix[j, k];
                            }
                }
                for (int i = 0; i < matrix.GetLength(0); i++)
                    if (check[i] < d[i]) { Console.WriteLine("Negative cicle"); return ; }
                
                Console.WriteLine("Non negative Cicle");
                Draw();
                Console.ReadLine();
                double psi = (2 * Math.PI) / matrix.GetLength(0);
                for (int i = 0; i < matrix.GetLength(0); i++)
                {

                    Picture.DrawNumbersInverse(d[i], R_bigCircle * Math.Cos(psi * i) + r_smallCircle, R_bigCircle * Math.Sin(psi * i) + r_smallCircle, 6);
                }
                Gl.glFlush();
                Console.ReadLine();


            }
            public bool _BellmanFord(int s,  int[] d)//or
            {
                //d = new int[matrix.GetLength(0)];
                for (int i = 0; i < d.Length; i++)
                    d[i] = INF;

                d[s] = 0;
                int l = 0;
                for (; l < matrix.GetLength(0); l++)
                {
                    for (int j = 0; j < matrix.GetLength(0); j++)
                        for (int k = 0; k < matrix.GetLength(1); k++)
                            if (matrix[j, k] != INF)
                            {
                                if (d[k] > d[j] + matrix[j, k]) d[k] = d[j] + matrix[j, k];
                            }
                }
                int[] check = new int[matrix.GetLength(0)];
                for (int i = 0; i < d.Length; i++)
                {
                    check[i] = d[i];
                }                for (; l < matrix.GetLength(0) + 1; l++)
                {
                    for (int j = 0; j < matrix.GetLength(0); j++)
                        for (int k = 0; k < matrix.GetLength(1); k++)
                            if (matrix[j, k] != INF)
                            {
                                if (check[k] > check[j] + matrix[j, k]) check[k] = check[j] + matrix[j, k];
                            }
                }
                for (int i = 0; i < matrix.GetLength(0); i++)
                    if (check[i] < d[i]) { Console.WriteLine("Negative cicle"); return false; }
                //Console.WriteLine("Non negative Cicle");
                return true;
            }
            public void _Deikstra(int num_v,  int [] d,  int[] p)
            {
                //d = new int[matrix.GetLength(0)];
               // p = new int[matrix.GetLength(0)];
                int[] sel_e = new int[matrix.GetLength(0)];
                bool[] used = new bool[matrix.GetLength(0)];

                for (int i = 0; i < matrix.GetLength(0); i++)
                {
                    d[i] = INF;
                    sel_e[i] = -1;
                    p[i] = -1;
                    used[i] = false;
                }


                d[num_v] = 0;

                for (int i = 0; i < matrix.GetLength(0); ++i)
                {
                    int v = -1;
                    for (int j = 0; j < matrix.GetLength(0); ++j)
                        if (!used[j] && (v == -1 || d[j] < d[v]))
                            v = j;
                    if (d[v] == INF) { Console.WriteLine("No MST!"); return; }

                    used[v] = true;
                    color_Vertex[v] = 3;

                    for (int k = 0; k < matrix.GetLength(0); ++k)
                    {
                        if ((matrix[v, k] != INF) && (!used[k]))
                        {
                            int to = k;
                            int len = matrix[v, k];
                            if (d[v] + len < d[to])
                            {
                                d[to] = d[v] + len;
                                p[to] = v;
                            }
                        }
                    }
                }
            
            }
            public void Johnson()//or
            {
               int [][] d=new int[matrix.GetLength(0)][]; 
               int[][] p=new int[matrix.GetLength(0)][];
               int[] h=new int[matrix.GetLength(0)+1] ;
               Graph _G = new Graph();
               Graph G = new Graph();

               _G.color_Vertex = new int[matrix.GetLength(0)+1];
               _G.color_Edge = new int[matrix.GetLength(0)+1, matrix.GetLength(1)+1];
               _G.matrix = new int[matrix.GetLength(0)+1, matrix.GetLength(1)+1];
               G.color_Vertex = new int[matrix.GetLength(0)];
               G.color_Edge = new int[matrix.GetLength(0), matrix.GetLength(1)];
               G.matrix = new int[matrix.GetLength(0), matrix.GetLength(1)];

               for (int i = 0; i <= matrix.GetLength(0); i++)
               {
                   _G.color_Vertex[i] = 1;


                   for (int j = 0; j <= matrix.GetLength(1); j++)
                   {
                       _G.color_Edge[i, j] = 2;
                       if (i < matrix.GetLength(0) && j < matrix.GetLength(1))
                           _G.matrix[i, j] = matrix[i, j];
                       else
                       {
                           if (i == matrix.GetLength(0)) _G.matrix[i, j] = 0;
                           else _G.matrix[i, j] = INF;
                       } //!!
                   }
               }
               for (int i = 0; i < matrix.GetLength(0); i++)
               {
                   d[i] = new int[matrix.GetLength(1)];
                   p[i] = new int[matrix.GetLength(1)];
                   G.color_Vertex[i] = 1;
                   for (int j = 0; j < matrix.GetLength(1); j++)
                   {
                       G.color_Edge[i, j] = 2;
                       G.matrix[i, j] = matrix[i, j];
                       
                   }
               }


               bool f = _G._BellmanFord(matrix.GetLength(0), h);
               if (f == false) return;// граф содержит цикл с отрицательным весом

             for (int i = 0; i < matrix.GetLength(0); i++)
                 for (int j = 0; j < matrix.GetLength(1); j++)
                 if( G.matrix[i,j]!=INF)
                 {
                     G.matrix[i, j] = matrix[i, j] + h[i] - h[j];
                   
                 }
        
             for (int k=0; k < matrix.GetLength(0); k++)
                 G._Deikstra(k,  d[k],  p[k]);
             for (int i = 0; i < matrix.GetLength(0); i++)
                 for (int j = 0; j < matrix.GetLength(1); j++)
                     if (G.matrix[i, j] != INF)
                     {
                         d[i][ j] = d[i][ j] - h[i] + h[j];
                         
                     }
             Console.WriteLine("Johnson algorithm\n\n");
             for (int k = 0; k < matrix.GetLength(0); ++k)
             {
                 for (int j = 0; j < matrix.GetLength(1); ++j)
                 { Console.Write(d[k][ j] + " "); }
                 Console.WriteLine();
             }

             Console.WriteLine();
            }

            public bool _BFS(int start, int end, int [,] capacity, int[,] flow, out int[] pred)
            {
                bool[] used = new bool[matrix.GetLength(0)];
                for (int i = 0; i < matrix.GetLength(0); i++)
                { used[i]= false; }
                used[start] = true;
                int t = 0, h = 0;
                int [] queue = new int[matrix.GetLength(0)];
                pred = new int[matrix.GetLength(0)];
                queue[t++] = start;
                pred[start] = -1;
                while (h < t)
                {
                    int v = queue[h]; h++;
                    for (int j = 0; j < matrix.GetLength(0); ++j)
                    {
                        if (matrix[v, j] != INF)
                        {
                            if ((!used[j])&&(capacity[v,j]-flow[v,j]>0))
                            {
                                used[j] = true;
                                queue[t++] = j;
                                pred[j] = v;

                            }
                        }
                    }
                }
                if (used[end]) return false;
                else return true;
            }
            public int FordFalkerson(int source, int stock)
            {
                int [,] capacity=new int [matrix.GetLength(0),matrix.GetLength(1)];//матрица пропускных способностей
                int[,] flow = new int[matrix.GetLength(0), matrix.GetLength(1)];//матрица потока
                //int [] pred = new int [matrix.GetLength(0)];// массив пути
                int [] pred;
                int maxflow = 0;
                for (int i = 0; i < matrix.GetLength(0); i++)
                    for (int j = 0; j < matrix.GetLength(1); j++)
                    { 
                        if(matrix[i,j]!=INF)
                        capacity[i,j]=matrix[i,j];
                        else capacity[i,j]=0;
                        flow[i, j] = 0;
                    }

                
                while (!_BFS(source,stock, capacity, flow, out pred))
                {
                    int delta = 10000;
                    for (int u = matrix.GetLength(0) - 1; pred[u] >= 0; u = pred[u])


                        delta = min(delta, (capacity[pred[u], u] - flow[pred[u], u]));

                    for (int u = matrix.GetLength(0) - 1; pred[u] >= 0; u = pred[u])
                    {
                        flow[pred[u], u] += delta;
                        flow[u, pred[u]] -= delta;
                        color_Edge[pred[u], u] = 4;
                        color_Edge[ u,pred[u]] = 4;
                    }
                    maxflow += delta;
                    Draw_flow(capacity, flow);
                    Console.ReadLine();
                    for (int i = 0; i < matrix.GetLength(0); i++)
                        for (int j = 0; j < matrix.GetLength(1); j++)
                        {
                            if ((capacity[i, j] != 0) && (capacity[i, j] == flow[i, j] || capacity[i, j] == -flow[i, j]))
                            {
                                color_Edge[i, j] = 6;
                                color_Edge[j , i] = 6;
                            }
                        }
                    Draw_flow(capacity, flow);
                }

                Picture.DrawNumbersInverse(maxflow, 0.8 , -0.8 , 2);
                Gl.glFlush();
                Console.ReadLine();
                return maxflow;
            }

            public bool _Deikstra(int num_v, int end, out int[] p, int [,] capacity, int[,] flow)
            {
                int [] d = new int[matrix.GetLength(0)];
                p = new int[matrix.GetLength(0)];
                
                int[] sel_e = new int[matrix.GetLength(0)];
                bool[] used = new bool[matrix.GetLength(0)];

                for (int i = 0; i < matrix.GetLength(0); i++)
                {
                    d[i] = INF;
                    sel_e[i] = -1;
                    p[i] = -1;
                    used[i] = false;
                }


                d[num_v] = 0;

                for (int i = 0; i < matrix.GetLength(0); ++i)
                {
                    int v = -1;
                    for (int j = 0; j < matrix.GetLength(0); ++j)
                        if (!used[j] && (v == -1 || d[j] < d[v]))
                            v = j;
                    if (d[v] == INF) { continue ; }

                    used[v] = true;
                   

                    for (int k = 0; k < matrix.GetLength(0); ++k)
                    {
                        if ((matrix[v, k] != INF) && (!used[k])&&(capacity[v,k]-flow[v,k]>0))
                        {
                            int to = k;
                            int len = matrix[v, k];
                            if (d[v] + len < d[to])
                            {
                                d[to] = d[v] + len;
                                p[to] = v;
                            }
                        }
                    }
                }
                if (used[end]) return false;
                else return true;

            }
            public int EdmondsKarp(int source, int stock)
            {

                int[,] capacity = new int[matrix.GetLength(0), matrix.GetLength(1)];//матрица пропускных способностей
                int[,] flow = new int[matrix.GetLength(0), matrix.GetLength(1)];//матрица потока
               
                int[] pred;//массив пути
                int maxflow = 0;
                for (int i = 0; i < matrix.GetLength(0); i++)
                    for (int j = 0; j < matrix.GetLength(1); j++)
                    {
                        if (matrix[i, j] != INF)
                            capacity[i, j] = matrix[i, j];
                        else capacity[i, j] = 0;
                        flow[i, j] = 0;
                    }


                while (!_Deikstra(source, stock, out pred, capacity, flow))
                {
                    int delta = 10000;
                    for (int u = matrix.GetLength(0) - 1; pred[u] >= 0; u = pred[u])


                        delta = min(delta, (capacity[pred[u], u] - flow[pred[u], u]));

                    for (int u = matrix.GetLength(0) - 1; pred[u] >= 0; u = pred[u])
                    {
                        flow[pred[u], u] += delta;
                        flow[u, pred[u]] -= delta;
                        color_Edge[pred[u], u] = 4;
                        color_Edge[u, pred[u]] = 4;

                    }
                    maxflow += delta;
                    Draw_flow(capacity, flow);
                    Console.ReadLine();
                    for (int i = 0; i < matrix.GetLength(0); i++)
                        for (int j = 0; j < matrix.GetLength(1); j++)
                        {
                            if ((capacity[i, j] != 0) && (capacity[i, j] == flow[i, j] || capacity[i, j] == -flow[i, j]))
                            {
                                color_Edge[i, j] = 6;
                                color_Edge[j, i] = 6;
                            }
                        }
                    Draw_flow(capacity, flow);
                }
                Picture.DrawNumbersInverse(maxflow, 0.8, -0.8, 2);
                Gl.glFlush();
                Console.ReadLine();
                return maxflow;
            }


            public void CreateListEdge(out List<edge> wEdge)
            {
                wEdge = new List<edge>();


                for (int i = 0; i < matrix.GetLength(0); i++)
                    for (int j = i; j < matrix.GetLength(1); j++) // для несимметрических j=0
                    {
                        if (matrix[i, j] != INF)
                        {
                            edge w = new edge();
                            w.from = i;
                            w.to = j;
                            w.weight = matrix[i, j];
                           // w.color = color_Edge[i, j];
                            wEdge.Add(w);
                        }
                    }
            }
            public void CreateListVertex(out List<vertex> vVertex)
            {
                vVertex = new List<vertex>();

                for (int i = 0; i < matrix.GetLength(0); i++)
                {
                    vertex v = new vertex();
                    for (int j = 0; j < matrix.GetLength(1); j++)
                    {

                        if (matrix[i, j] != INF)
                        {
                            edge w = new edge();
                            w.from = i;
                            w.to = j;
                            w.weight = matrix[i, j];
                           // w.color = color_Edge[i, j];
                            v.rebra.Add(w);
                        }
                    }
                    v.number = i;
                    v.color = color_Vertex[i];
                    vVertex.Add(v);
                }
            }
            public void SortByWeightListEdge(List<edge> wEdge)
            {
                wEdge.Sort(delegate(edge w1, edge w2)
                { return w1.weight.CompareTo(w2.weight); });
            }
            public void FromListOfEdgesToGraph(List<edge> wEdge)
            {
               // Console.WriteLine("\nFromListOfEdgesToGraph\n");
                FromListOfEdgesToMatrix(wEdge, out matrix);
                DeleteZeroLineNColumn(ref matrix);

                color_Edge = new int[matrix.GetLength(0), matrix.GetLength(0)];
                color_Vertex = new int[matrix.GetLength(0)];
                for (int i = 0; i < matrix.GetLength(0); i++) color_Vertex[i] = 1;
                for (int j = 0; j < matrix.GetLength(0); j++)
                    for (int k = 0; k < matrix.GetLength(1); k++)
                    {
                        if (matrix[j, k] != 0) color_Edge[j, k] = 2;
                        else matrix[j, k] = INF;
                    }

                //Console.Write(" \nVertex\n");
                //foreach (var element in color_Vertex) Console.Write(element + " ");
                //Console.Write(" \nEdge\n");

                //for (int k = 0; k < color_Edge.GetLength(0); ++k)
                //{
                //    for (int j = 0; j < color_Edge.GetLength(1); ++j)
                //    { Console.Write(color_Edge[k, j] + " "); }
                //    Console.WriteLine();
                //}

                Console.WriteLine();
            }// не переносит цвета
            public void FromListOfVertexsToGraph(List<vertex> vVertex)
            {
               // Console.WriteLine("\nFromListOfVertexsToGraph\n");
                List<edge> wEdge;
                FromListOfVertexsToListOfEdges(vVertex, out wEdge);
                FromListOfEdgesToMatrix(wEdge, out matrix);
                DeleteZeroLineNColumn(ref matrix);

                color_Edge = new int[matrix.GetLength(0), matrix.GetLength(0)];
                color_Vertex = new int[matrix.GetLength(0)];
                for (int i = 0; i < matrix.GetLength(0); i++) color_Vertex[i] = 1;
                for (int j = 0; j < matrix.GetLength(0); j++)
                    for (int k = 0; k < matrix.GetLength(1); k++)
                    {
                        if (matrix[j, k] != 0) color_Edge[j, k] = 2;
                        else matrix[j, k] = INF;
                    }

                Console.Write(" \nVertex\n");
                foreach (var element in color_Vertex) Console.Write(element + " ");
                Console.Write(" \nEdge\n");

                for (int k = 0; k < color_Edge.GetLength(0); ++k)
                {
                    for (int j = 0; j < color_Edge.GetLength(1); ++j)
                    { Console.Write(color_Edge[k, j] + " "); }
                    Console.WriteLine();
                }

                Console.WriteLine();
            }// не переносит цвета


        }

        static void FromListOfEdgesToMatrix(List<edge> wEdge, out int[,] massive)
        {
            massive = new int[NMax, NMax];
            for (int i = 0; i < wEdge.Count; i++)
            {
                massive[wEdge[i].from, wEdge[i].to] = wEdge[i].weight;
                massive[wEdge[i].to, wEdge[i].from] = wEdge[i].weight;
            }
            //Console.WriteLine(" \nMassive:\n");
            //for (int k = 0; k < massive.GetLength(0); ++k)
            //{
            //    for (int j = 0; j < massive.GetLength(1); ++j)
            //    { Console.Write(massive[k, j] + " "); }
            //    Console.WriteLine();
            //}
        }
        static void FromListOfVertexsToListOfEdges(List<vertex> vVertex, out List<edge> wEdge)
        {
            wEdge = new List<edge>();

            for (int i = 0; i < vVertex.Count; i++)
            {
                wEdge.AddRange(vVertex[i].rebra);
            }

            for (int j = 0; j < wEdge.Count; j++)
            {
                edge w = new edge();
                w = wEdge[j];
                for (int k = 0; k < wEdge.Count; k++)
                {
                    if ((w.weight == wEdge[k].weight) && (w.from == wEdge[k].to) && (w.to == wEdge[k].from))
                    {
                        wEdge.RemoveAt(k);
                    }
                }
            }
        }
        static int min(int a, int b)
        {
            if (a < b) return a;
            else return b;
        }
        static void DeleteZeroLineNColumn(ref int[,] massive)
        {
            int[,] temp = new int[massive.GetLength(0), massive.GetLength(1)];
            int rank = 0;
            for (int i = 0; i < massive.GetLength(0); i++)
            {
                int k = 0;
                for (int j = 0; j < massive.GetLength(1); j++)
                {
                    if (massive[i, j] == 0) k++;
                }
                if (k != massive.GetLength(1))
                {
                    for (int m = 0; m < massive.GetLength(1); m++)
                    {
                        temp[rank, m] = massive[i, m];
                    }
                    rank++;
                }
            }
            int[,] res = new int[rank, rank];
            int n = 0;
            for (int j = 0; j < temp.GetLength(1); j++)
            {
                int k = 0;
                for (int i = 0; i < rank; i++)
                {
                    if (temp[i, j] == 0) k++;
                }
                if (k != rank)
                {
                    for (int m = 0; m < rank; m++)
                    {
                        res[m, n] = temp[m, j];
                    }
                    n++;
                }
            }
            massive = res;
            //Console.WriteLine(" \n Result:\n");
            //for (int i = 0; i < massive.GetLength(0); ++i)
            //{
            //    for (int j = 0; j < massive.GetLength(1); ++j)
            //    { Console.Write(massive[i, j] + " "); }
            //    Console.WriteLine();
            //}

        }
        static void ReadFromFileLineByLineIntoArray(ref int[,] arrMain2)
        {
            string inputString;

            StreamReader file = new StreamReader(@"C:\Users\user01\Documents\Visual Studio 2012\Projects\Graph_Console\do.txt");//!!
            int i = 0;
            string temp = "";

            int[,] arrMain = new int[NMax, NMax];
            while ((inputString = file.ReadLine()) != null)
            {
                string[] spl = null;
                int j = 0;
                spl = inputString.Split(new Char[] { ' ', '|' });

                foreach (string Mas in spl)
                    if (!(temp.Equals(Mas)))
                        arrMain[i, j++] = Convert.ToInt32(Mas);
                ++i;
            }
            file.Close();

            int[,] arrMainLast = new int[i, i];

            for (int k = 0; k < arrMainLast.GetLength(0); ++k)
            {
                for (int j = 0; j < arrMainLast.GetLength(1); ++j)
                {
                    arrMainLast[k, j] = arrMain[k, j];
                    Console.Write(arrMainLast[k, j] + " ");
                }
                Console.WriteLine();
            }
            arrMain2 = arrMainLast;
        }
        static void FromEdgeToMatrix(edge w, out int[,] massive)
        {
            massive = new int[NMax, NMax];
            for (int i = 0; i < NMax; i++)
            {
                massive[w.from, w.to] = w.weight;
                massive[w.to, w.from] = w.weight;
            }

        }

        public class Picture
        {
            // функция визуализации текста 
            public static void PrintText2D(float x, float y, string text, string str = "blue")
            {
                if (str == "red")
                    Gl.glColor3f(1, 0, 0);
                if (str == "green")
                    Gl.glColor3f(0, 1, 0);
                if (str == "blue")
                    Gl.glColor3f(0, 0, 1);
                // устанавливаем позицию вывода растровых символов 
                // в переданных координатах x и y. 
                Gl.glRasterPos2f(x, y);

                // в цикле foreach перебираем значения из массива text, 
                // который содержит значение строки для визуализации 
                foreach (char char_for_draw in text)
                {
                    // символ C визуализируем с помощью функции glutBitmapCharacter, используя шрифт GLUT_BITMAP_9_BY_15. 
                    Glut.glutBitmapCharacter(Glut.GLUT_BITMAP_9_BY_15, char_for_draw);
                }

            }
            public static void Line(double x1, double y1, double x2, double y2, int n)
            {
                fColor(n);

                Gl.glBegin(Gl.GL_LINE_STRIP);

                Gl.glVertex2d(x1, y1);
                Gl.glVertex2d(x2, y2);

                Gl.glEnd();


            }
            public static void Circle(double x, double y, int n)
            {


                fColor(n);

                const double h = 0.005;
                Gl.glBegin(Gl.GL_TRIANGLES);

                double i = 0.0;
                while (i <= 2 * Math.PI)
                {
                    Gl.glVertex2d(x, y);
                    Gl.glVertex2d(x + r_smallCircle * Math.Cos(i), y + r_smallCircle * Math.Sin(i));
                    Gl.glVertex2d(x + r_smallCircle * Math.Cos(i + h), y + r_smallCircle * Math.Sin(i + h));
                    i = i + h;
                }
                Gl.glEnd();


            }


            public static void fColor(int n)
            {

                switch (n)
                {
                    case 0:
                        {
                            Gl.glColor3f(0, 0, 0);//black
                            break;
                        }
                    case 1:
                        {
                            Gl.glColor3f(255, 255, 0);//yellow
                            break;
                        }
                    case 2:
                        {
                            Gl.glColor3f(255, 0, 0);//red
                            break;
                        }
                    case 3:
                        {
                            Gl.glColor3f(0, 0, 255);//blue
                            break;
                        }
                    case 4:
                        {
                            Gl.glColor3f(0, 255, 0);//green
                            break;
                        }
                    case 5:
                        {

                            Gl.glColor3f(0, 255, 255);//голубой
                            break;
                        }
                    case 6:
                        {
  
                            Gl.glColor3f(255, 0, 255);//пурпурный
                            break;
                        }
                    default:
                        {
                            Gl.glColor3f(0, 0, 0);//black
                            break;
                        }

                }
            }
            public static void Numerals(int numeral, double x, double y)
            {

                switch (numeral)
                {
                    case 0:
                        {
                            Gl.glBegin(Gl.GL_LINE_LOOP);
                            Gl.glVertex2d(x - h_num / 2, y + h_num);
                            Gl.glVertex2d(x + h_num / 2, y + h_num);
                            Gl.glVertex2d(x + h_num / 2, y - h_num);
                            Gl.glVertex2d(x - h_num / 2, y - h_num);
                            Gl.glEnd();
                            break;
                        }
                    case 1:
                        {
                            Gl.glBegin(Gl.GL_LINE_STRIP);

                            Gl.glVertex2d(x + h_num / 2, y + h_num);
                            Gl.glVertex2d(x + h_num / 2, y - h_num);

                            Gl.glEnd();
                            break;
                        }
                    case 2:
                        {
                            Gl.glBegin(Gl.GL_LINE_STRIP);

                            Gl.glVertex2d(x - h_num / 2, y + h_num);
                            Gl.glVertex2d(x + h_num / 2, y + h_num);
                            Gl.glVertex2d(x + h_num / 2, y);
                            Gl.glVertex2d(x - h_num / 2, y);
                            Gl.glVertex2d(x - h_num / 2, y - h_num);
                            Gl.glVertex2d(x + h_num / 2, y - h_num);

                            Gl.glEnd();
                            break;
                        }
                    case 3:
                        {
                            Gl.glBegin(Gl.GL_LINE_STRIP);

                            Gl.glVertex2d(x - h_num / 2, y + h_num);
                            Gl.glVertex2d(x + h_num / 2, y + h_num);
                            Gl.glVertex2d(x + h_num / 2, y - h_num);
                            Gl.glVertex2d(x - h_num / 2, y - h_num);

                            Gl.glEnd();

                            Gl.glBegin(Gl.GL_LINES);

                            Gl.glVertex2d(x - h_num / 2, y);
                            Gl.glVertex2d(x + h_num / 2, y);

                            Gl.glEnd();

                            break;
                        }
                    case 4:
                        {
                            Gl.glBegin(Gl.GL_LINE_STRIP);

                            Gl.glVertex2d(x - h_num / 2, y + h_num);
                            Gl.glVertex2d(x - h_num / 2, y);
                            Gl.glVertex2d(x + h_num / 2, y);
                            Gl.glVertex2d(x + h_num / 2, y - h_num);

                            Gl.glEnd();

                            Gl.glBegin(Gl.GL_LINES);

                            Gl.glVertex2d(x + h_num / 2, y + h_num);
                            Gl.glVertex2d(x + h_num / 2, y);

                            Gl.glEnd();

                            break;
                        }
                    case 5:
                        {
                            Gl.glBegin(Gl.GL_LINE_STRIP);

                            Gl.glVertex2d(x + h_num / 2, y + h_num);
                            Gl.glVertex2d(x - h_num / 2, y + h_num);
                            Gl.glVertex2d(x - h_num / 2, y);
                            Gl.glVertex2d(x + h_num / 2, y);
                            Gl.glVertex2d(x + h_num / 2, y - h_num);
                            Gl.glVertex2d(x - h_num / 2, y - h_num);

                            Gl.glEnd();
                            break;
                        }
                    case 6:
                        {
                            Gl.glBegin(Gl.GL_LINE_STRIP);

                            Gl.glVertex2d(x + h_num / 2, y + h_num);
                            Gl.glVertex2d(x - h_num / 2, y + h_num);
                            Gl.glVertex2d(x - h_num / 2, y);
                            Gl.glVertex2d(x + h_num / 2, y);
                            Gl.glVertex2d(x + h_num / 2, y - h_num);
                            Gl.glVertex2d(x - h_num / 2, y - h_num);
                            Gl.glVertex2d(x - h_num / 2, y);

                            Gl.glEnd();
                            break;
                        }
                    case 7:
                        {
                            Gl.glBegin(Gl.GL_LINE_STRIP);

                            Gl.glVertex2d(x - h_num / 2, y + h_num);
                            Gl.glVertex2d(x + h_num / 2, y + h_num);
                            Gl.glVertex2d(x + h_num / 2, y - h_num);

                            Gl.glEnd();
                            break;
                        }
                    case 8:
                        {
                            Gl.glBegin(Gl.GL_LINE_STRIP);

                            Gl.glVertex2d(x + h_num / 2, y);
                            Gl.glVertex2d(x - h_num / 2, y);
                            Gl.glVertex2d(x - h_num / 2, y + h_num);
                            Gl.glVertex2d(x + h_num / 2, y + h_num);
                            Gl.glVertex2d(x + h_num / 2, y - h_num);
                            Gl.glVertex2d(x - h_num / 2, y - h_num);
                            Gl.glVertex2d(x - h_num / 2, y);
                            Gl.glEnd();

                            break;
                        }
                    case 9:
                        {
                            Gl.glBegin(Gl.GL_LINE_STRIP);

                            Gl.glVertex2d(x - h_num / 2, y - h_num);
                            Gl.glVertex2d(x + h_num / 2, y - h_num);
                            Gl.glVertex2d(x + h_num / 2, y + h_num);
                            Gl.glVertex2d(x - h_num / 2, y + h_num);
                            Gl.glVertex2d(x - h_num / 2, y);
                            Gl.glVertex2d(x + h_num / 2, y);

                            Gl.glEnd();

                            break;
                        }
                    default:
                        {
                            break;
                        }

                }
            }
            public static void DrawNumerals(int numeral, double x, double y)
            {

                Gl.glClear(Gl.GL_COLOR_BUFFER_BIT | Gl.GL_DEPTH_BUFFER_BIT);
                Picture.fColor(0);
                Numerals(numeral, x, y);
                Gl.glFlush();
            }
            public static void DrawNumbers(int number, double x, double y)
            {
                Picture.fColor(0);
                
                if (number == 0)
                { 
                    Numerals(number, x, y);
             
                    return;
                }

                int n = number, k = 0;
                const double step = 0.01;


                while (n > 0)
                {
                    n = n / 10; ++k;
                }
                n = number;
                if (k % 2 == 0)
                {
                    for (int i = 0; n > 0; i++)
                    {
                        Numerals(n % 10, x - ((k - 1) * h_num) / 2 + i * (h_num + step), y);
                        n = n / 10;
                    }
                }
                if (k % 2 != 0)
                {
                    for (int i = 0; n > 0; i++)
                    {
                        Numerals(n % 10, x - (k * h_num) / 2 + i * (h_num + step), y);
                        n = n / 10;
                    }
                }

            }
            public static void DrawNumbersInverse(int number, double x, double y, int color)
            {
                Picture.fColor(color);
                if (number == 0)
                { Numerals(number, x, y); return; }

                int n = number, k = 0;
                const double step = 0.01;

                while (n > 0)
                {
                    n = n / 10; ++k;
                }
                n = Math.Abs(number);
                if (number > 0)
                {
                    if (k % 2 == 0)
                    {
                        for (int i = 0; n > 0; i++)
                        {
                            Numerals(n % 10, x + ((k - 1) * h_num) / 2 - i * (h_num + step), y);
                            n = n / 10;
                        }
                    }
                    if (k % 2 != 0)
                    {
                        for (int i = 0; n > 0; i++)
                        {
                            Numerals(n % 10, x + (k * h_num) / 2 - i * (h_num + step), y);
                            n = n / 10;
                        }
                    }

                }

                else
                {
                    if (k % 2 == 0)
                    {
                        for (int i = 0; n > 0; i++)
                        {
                            if (i == 0)
                            {
                                Gl.glBegin(Gl.GL_LINES);

                                Gl.glVertex2d(x + ((k-1) * h_num) / 2 - i * (h_num + step) - (3*h_num) / 2-0.01, y);
                                Gl.glVertex2d(x + ((k-1) * h_num) / 2 - i * (h_num + step)-0.01, y);

                                Gl.glEnd();
                            }
                            Numerals(n % 10, x + ((k-1) * h_num) / 2 - i * (h_num + step), y);
                            n = n / 10;
                        }
                    }
                    if (k % 2 != 0)
                    {

                        for (int i = 0; n > 0; i++)
                        {
                            if (i == 0)
                            {
                                Gl.glBegin(Gl.GL_LINES);

                                Gl.glVertex2d(x + ((k) * h_num) / 2 - i * (h_num + step) - (3*h_num / 2)-0.01, y);
                                Gl.glVertex2d(x + ((k) * h_num) / 2 - i * (h_num + step) -0.01, y);

                                Gl.glEnd();
                            }
                            Numerals(n % 10, x + ((k) * h_num) / 2 - i * (h_num + step), y);
                            n = n / 10;
                        }
                    }

                }
            }


        }


        static void Display()
        {
           A.Draw();               
        }
        static void Update()
        {
            Glut.glutPostRedisplay();
        }

        static void KeyboardEvent(int key, int x, int y)
        {
            switch (key)
            {
                case Glut.GLUT_KEY_F1:
                    {
                        A.First();
                        A.Draw();
                        Console.ReadLine();
                        string text="DFS";
                        Console.WriteLine(text);
                        //Picture.PrintText2D(0.0f,0.0f,text,"blue");
                        A.DFS(0);
                        Console.ReadLine();


                        break;
                    }
                case Glut.GLUT_KEY_F2:
                    {
                        A.First();
                        A.Draw();
                        Console.ReadLine();
                        Console.WriteLine("BFS");
                        A.BFS(0);
                        Console.ReadLine();


                        break;
                    }
                case Glut.GLUT_KEY_F3:
                    {
                        A.First();
                        A.Draw();
                        Console.ReadLine();
                        Console.WriteLine("Kruskala");
                        A.Kruskala(); Console.ReadLine();

                        break;
                    }
                case Glut.GLUT_KEY_F4:
                    {
                        A.First();
                        A.Draw();
                        Console.ReadLine();
                        Console.WriteLine("Prima");
                        A.Prima(0); Console.ReadLine();

                        break;
                    }
                case Glut.GLUT_KEY_F5:
                    {
                        A.First();
                        A.Draw();
                        Console.ReadLine();
                        Console.WriteLine("BellmanFord");
                        A.BellmanFord(0); Console.ReadLine();

                        break;
                    }
                case Glut.GLUT_KEY_F6:
                    {
                        A.First();
                        A.Draw();
                        Console.ReadLine();
                        Console.WriteLine("Deikstra");
                        A.Deikstra(0); Console.ReadLine();

                        break;
                    }
                case Glut.GLUT_KEY_F7:
                    {
                        A.First();
                        A.Draw();
                        Console.ReadLine();
                        Console.WriteLine("FloydWarshall");
                        A.FloydWarshall(); Console.ReadLine();

                        break;
                    }
                case Glut.GLUT_KEY_F8:
                    {
                        A.First();
                        A.Draw();
                        Console.ReadLine();
                        Console.WriteLine("Johnson");
                        A.Johnson(); Console.ReadLine();

                        break;
                    }
                case Glut.GLUT_KEY_F9:
                    {
                        A.First();
                        A.Draw();
                        Console.ReadLine();
                        Console.WriteLine("FordFalkerson");
                        A.FordFalkerson(0, A.matrix.GetLength(0) - 1); Console.ReadLine();

                        break;
                    }
                case Glut.GLUT_KEY_F10:
                    {
                        A.First();
                        A.Draw();
                        Console.ReadLine();
                        Console.WriteLine("EdmondsKarp");
                        A.EdmondsKarp(0, A.matrix.GetLength(0) - 1); Console.ReadLine();

                        break;
                    }
            }
        }


        static void Main(string[] args)
        {
            A.Initialize();

            Glut.glutInit();
            Glut.glutInitWindowSize(600, 600);
            Glut.glutInitWindowPosition(0, 0);
            Glut.glutCreateWindow(" Tao Graph");

            Gl.glClearColor(255, 255, 255, 1);

            //Glut.glutDisplayFunc(
            //   Display
            //);
           
            //Glut.glutIdleFunc(
            //        Update
            //    );
            
            Glut.glutSpecialFunc(
                             KeyboardEvent
            );






            Glut.glutMainLoop();
        }
    }
}








