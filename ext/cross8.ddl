* cdd+: Double Description Method in C++:Version 0.73 (September 6, 1995)
* Copyright (C) 1995, Komei Fukuda, fukuda@ifor.math.ethz.ch
* Compiled for Floating-Point Arithmetic
*Input File: cross8.ine(256 x 9)
*Initial set of hyperplanes:  64 80 84 85 86 87 89 97 129
begin
 247  5
10 193 11 10 8
11 161 12 10 8
12 32 15 12 8
13 225 15 11 8
14 145 16 11 8
15 48 18 12 8
16 209 18 11 8
17 113 22 14 8
18 177 22 13 8
19 16 23 13 8
20 241 23 12 8
21 137 24 12 8
22 56 26 13 8
23 201 26 12 8
24 105 29 14 8
25 169 29 13 8
26 24 30 13 8
27 233 30 12 8
28 72 35 16 8
29 153 35 15 8
30 40 36 15 8
31 217 36 14 8
32 121 38 15 8
33 185 38 14 8
34 8 39 14 8
35 249 39 13 8
36 133 40 13 8
37 60 42 14 8
38 197 42 13 8
39 101 45 15 8
40 165 45 14 8
41 28 46 14 8
42 229 46 13 8
43 76 50 16 8
44 149 50 15 8
45 44 51 15 8
46 213 51 14 8
47 117 53 15 8
48 181 53 14 8
49 12 54 14 8
50 245 54 13 8
51 93 60 18 8
52 141 60 17 8
53 52 61 17 8
54 205 61 16 8
55 109 63 17 8
56 173 63 16 8
57 20 64 16 8
58 237 64 15 8
59 68 67 17 8
60 157 67 16 8
61 36 68 16 8
62 221 68 15 8
63 125 70 16 8
64 189 70 15 8
65 4 71 15 8
66 253 71 14 8
67 131 72 14 8
68 62 74 15 8
69 195 74 14 8
70 99 77 16 8
71 163 77 15 8
72 30 78 15 8
73 227 78 14 8
74 78 82 17 8
75 147 82 16 8
76 46 83 16 8
77 211 83 15 8
78 115 85 16 8
79 179 85 15 8
80 14 86 15 8
81 243 86 14 8
82 91 91 18 8
83 139 91 17 8
84 54 92 17 8
85 203 92 16 8
86 107 94 17 8
87 171 94 16 8
88 22 95 16 8
89 235 95 15 8
90 70 98 17 8
91 155 98 16 8
92 38 99 16 8
93 219 99 15 8
94 123 101 16 8
95 187 101 15 8
96 6 102 15 8
97 251 102 14 8
98 82 109 20 8
99 135 109 19 8
100 58 110 19 8
101 199 110 18 8
102 103 112 19 8
103 167 112 18 8
104 26 113 18 8
105 231 113 17 8
106 74 116 19 8
107 151 116 18 8
108 42 117 18 8
109 215 117 17 8
110 119 119 18 8
111 183 119 17 8
112 10 120 17 8
113 247 120 16 8
114 95 124 19 8
115 143 124 18 8
116 50 125 18 8
117 207 125 17 8
118 111 127 18 8
119 175 127 17 8
120 18 128 17 8
121 239 128 16 8
122 66 131 18 8
123 159 131 17 8
124 34 132 17 8
125 223 132 16 8
126 127 134 17 8
127 191 134 16 8
128 2 135 16 8
129 255 135 15 8
130 130 136 15 9
131 63 138 16 10
132 194 138 15 10
133 98 141 17 11
134 162 141 16 11
135 31 142 16 11
136 226 142 15 11
137 79 146 18 12
138 146 146 17 12
139 47 147 17 12
140 210 147 16 12
141 114 149 17 12
142 178 149 16 12
143 15 150 16 12
144 242 150 15 12
145 90 155 19 13
146 138 155 18 13
147 55 156 18 13
148 202 156 17 13
149 106 158 18 13
150 170 158 17 13
151 23 159 17 13
152 234 159 16 13
153 71 162 18 13
154 154 162 17 13
155 39 163 17 13
156 218 163 16 13
157 122 165 17 13
158 186 165 16 13
159 7 166 16 13
160 250 166 15 13
161 83 172 20 14
162 134 172 19 14
163 59 173 19 14
164 198 173 18 14
165 102 175 19 14
166 166 175 18 14
167 27 176 18 14
168 230 176 17 14
169 75 179 19 14
170 150 179 18 14
171 43 180 18 14
172 214 180 17 14
173 118 182 18 14
174 182 182 17 14
175 11 183 17 14
176 246 183 16 14
177 94 187 19 14
178 142 187 18 14
179 51 188 18 14
180 206 188 17 14
181 110 190 18 14
182 174 190 17 14
183 19 191 17 14
184 238 191 16 14
185 67 194 18 14
186 158 194 17 14
187 35 195 17 14
188 222 195 16 14
189 126 197 17 14
190 190 197 16 14
191 3 198 16 14
192 254 198 15 14
193 88 206 22 16
194 132 206 21 16
195 61 207 21 16
196 196 207 20 16
197 100 209 21 16
198 164 209 20 16
199 29 210 20 16
200 228 210 19 16
201 77 213 21 16
202 148 213 20 16
203 45 214 20 16
204 212 214 19 16
205 116 216 20 16
206 180 216 19 16
207 13 217 19 16
208 244 217 18 16
209 92 221 21 16
210 140 221 20 16
211 53 222 20 16
212 204 222 19 16
213 108 224 20 16
214 172 224 19 16
215 21 225 19 16
216 236 225 18 16
217 69 228 20 16
218 156 228 19 16
219 37 229 19 16
220 220 229 18 16
221 124 231 19 16
222 188 231 18 16
223 5 232 18 16
224 252 232 17 16
225 81 237 21 16
226 136 237 20 16
227 57 238 20 16
228 200 238 19 16
229 104 240 20 16
230 168 240 19 16
231 25 241 19 16
232 232 241 18 16
233 73 244 20 16
234 152 244 19 16
235 41 245 19 16
236 216 245 18 16
237 120 247 19 16
238 184 247 18 16
239 9 248 18 16
240 248 248 17 16
241 96 252 20 16
242 144 252 19 16
243 49 253 19 16
244 208 253 18 16
245 112 255 19 16
246 176 255 18 16
247 17 256 18 16
248 240 256 17 16
249 65 259 19 16
250 160 259 18 16
251 33 260 18 16
252 224 260 17 16
253 128 262 18 16
254 192 262 17 16
255 1 263 17 16
256 256 263 16 16
end
*HyperplaneOrder: LexMin
*Degeneracy preknowledge for computation: None (possible degeneracy)
*Hull computation is chosen.
*set_intersection total#, effective#, loss# = 2150  488  1662
*Computation starts     at Fri Sep  8 17:16:02 1995
*            terminates at Fri Sep  8 17:16:11 1995
*Total processor time = 9 seconds
*                     = 0h 0m 9s
