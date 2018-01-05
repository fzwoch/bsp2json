/*
 * bsp2json
 *
 * Copyright (C) 2016 Florian Zwoch <fzwoch@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

package main

import (
	"encoding/binary"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"path/filepath"
	"strconv"
	"strings"

	"image"
	"image/color"
	"image/png"
)

type Entry struct {
	Offset uint32
	Size   uint32
}

type BspHeader struct {
	Version   uint32
	Entities  Entry
	Planes    Entry
	Miptex    Entry
	Vertices  Entry
	Visilist  Entry
	Nodes     Entry
	Texinfo   Entry
	Faces     Entry
	Lightmaps Entry
	Clipnodes Entry
	Leaves    Entry
	FacesList Entry
	Edges     Entry
	EdgesList Entry
	Models    Entry
}

type Vec3 struct {
	X float32
	Y float32
	Z float32
}

type BoundBox struct {
	Min Vec3
	Max Vec3
}

type Model struct {
	Bound     BoundBox
	Origin    Vec3
	NodeId0   uint32
	NodeId1   uint32
	NodeId2   uint32
	NodeId3   uint32
	Numleaves uint32
	FaceId    uint32
	FaceNum   uint32
}

type Miptex struct {
	Name    [16]byte
	Width   uint32
	Height  uint32
	Offset1 uint32
	Offset2 uint32
	Offset4 uint32
	Offset8 uint32
}

type Surface struct {
	VectorS   Vec3
	DistS     float32
	VectorT   Vec3
	DistT     float32
	TextureId uint32
	Animated  uint32
}

type Edge struct {
	Vertex0 uint16
	Vertex1 uint16
}

type Face struct {
	PlaneId   uint16
	Side      uint16
	LedgeId   int32
	LedgeNum  uint16
	TexinfoId uint16
	Typelight uint8
	Baselight uint8
	Light     [2]uint8
	Lightmap  int32
}

type UvStruct struct {
	u [3]float32
	v [3]float32
}

type FaceStruct struct {
	vert [3]int
	mat  int
}

var material_map map[uint32]int
var vertice_map map[uint16]int

var material_list []uint32
var vertice_list []uint16
var uv_list []UvStruct
var uv_light_list []int
var face_list []FaceStruct

func map_material(mat uint32) int {
	if i, ok := material_map[mat]; ok == true {
		return i
	}

	i := len(material_map)
	material_map[mat] = i
	material_list = append(material_list, mat)

	return i
}

func map_vertex(vert uint16) int {
	if i, ok := vertice_map[vert]; ok == true {
		return i
	}

	i := len(vertice_map)
	vertice_map[vert] = i
	vertice_list = append(vertice_list, vert)

	return i
}

func main() {
	if len(os.Args) < 3 {
		log.Fatalf("usage: %v <input.bsp> <output.json>", os.Args[0])
	}

	file, err := os.Open(os.Args[1])
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	var bsp BspHeader

	err = binary.Read(file, binary.LittleEndian, &bsp)
	if err != nil {
		log.Fatal(err)
	}

	if bsp.Version != 29 {
		log.Fatalf("BSP version %v != 29", bsp.Version)
	}

	var miptex_num uint32

	file.Seek(int64(bsp.Miptex.Offset), io.SeekStart)
	_ = binary.Read(file, binary.LittleEndian, &miptex_num)

	miptex_offsets := make([]uint32, miptex_num)
	miptex := make([]Miptex, miptex_num)

	_ = binary.Read(file, binary.LittleEndian, &miptex_offsets)

	for i := range miptex {
		file.Seek(int64(bsp.Miptex.Offset+miptex_offsets[i]), io.SeekStart)
		_ = binary.Read(file, binary.LittleEndian, &miptex[i])

		if miptex[i].Name[0] == '*' {
			miptex[i].Name[0] = '+'
		}

		if miptex[i].Name[0] == 0 {
			copy(miptex[i].Name[:], "unnamed"+strconv.Itoa(i))
		}
	}

	vertices := make([]Vec3, bsp.Vertices.Size/12)
	surfaces := make([]Surface, bsp.Texinfo.Size/40)
	faces := make([]Face, bsp.Faces.Size/20)
	lightmaps := make([]uint8, bsp.Lightmaps.Size)
	edges := make([]Edge, bsp.Edges.Size/4)
	edges_list := make([]int32, bsp.EdgesList.Size/4)
	models := make([]Model, bsp.Models.Size/64)

	_, _ = file.Seek(int64(bsp.Vertices.Offset), io.SeekStart)
	_ = binary.Read(file, binary.LittleEndian, &vertices)

	_, _ = file.Seek(int64(bsp.Texinfo.Offset), io.SeekStart)
	_ = binary.Read(file, binary.LittleEndian, &surfaces)

	_, _ = file.Seek(int64(bsp.Faces.Offset), io.SeekStart)
	_ = binary.Read(file, binary.LittleEndian, &faces)

	_, _ = file.Seek(int64(bsp.Lightmaps.Offset), io.SeekStart)
	_ = binary.Read(file, binary.LittleEndian, &lightmaps)

	_, _ = file.Seek(int64(bsp.Edges.Offset), io.SeekStart)
	_ = binary.Read(file, binary.LittleEndian, &edges)

	_, _ = file.Seek(int64(bsp.EdgesList.Offset), io.SeekStart)
	_ = binary.Read(file, binary.LittleEndian, &edges_list)

	_, _ = file.Seek(int64(bsp.Models.Offset), io.SeekStart)
	_ = binary.Read(file, binary.LittleEndian, &models)

/*
	lm_file, _ := os.Create("lightmap.png")

	var vv = int(math.Sqrt(float64(bsp.Lightmaps.Size)))
	lm_image := image.NewGray(image.Rect(0, 0, vv, vv))

	for i := 0; i < vv; i++ {
		for j := 0; j < vv; j++ {
			lm_image.Set(i, j, color.Gray{lightmaps[i*vv +j]})
		}
	}

	png.Encode(lm_file, lm_image)
	*/

	var lightmap_images []*image.Gray

	for model_idx, model := range models {
		material_map = make(map[uint32]int)
		vertice_map = make(map[uint16]int)

		material_list = nil
		vertice_list = nil
		uv_list = nil
		uv_light_list = nil
		face_list = nil

		for _, face := range faces[model.FaceId : model.FaceId+model.FaceNum] {
			var u_min float64 = math.MaxFloat64
			var u_max float64
			var v_min float64 = math.MaxFloat64
			var v_max float64

			var vertex = make([]uint16, face.LedgeNum)

			for i := range vertex {
				edge_id := edges_list[face.LedgeId+int32(i)]

				if edge_id < 0 {
					edge := edges[-edge_id]
					vertex[i] = edge.Vertex1
				} else {
					edge := edges[edge_id]
					vertex[i] = edge.Vertex0
				}
			}

			surface := surfaces[face.TexinfoId]
			mip := miptex[surface.TextureId]

			for i := 1; i < len(vertex)-1; i++ {
				var face FaceStruct
				var uv UvStruct

				face.vert[0] = map_vertex(vertex[0])
				face.vert[1] = map_vertex(vertex[i+1])
				face.vert[2] = map_vertex(vertex[i])
				face.mat = map_material(surface.TextureId)

				face_list = append(face_list, face)

				uv.u[0] = vertices[vertex[0]].X*surface.VectorS.X + vertices[vertex[0]].Y*surface.VectorS.Y + vertices[vertex[0]].Z*surface.VectorS.Z + surface.DistS
				uv.v[0] = vertices[vertex[0]].X*surface.VectorT.X + vertices[vertex[0]].Y*surface.VectorT.Y + vertices[vertex[0]].Z*surface.VectorT.Z + surface.DistT

				uv.u[2] = vertices[vertex[i]].X*surface.VectorS.X + vertices[vertex[i]].Y*surface.VectorS.Y + vertices[vertex[i]].Z*surface.VectorS.Z + surface.DistS
				uv.v[2] = vertices[vertex[i]].X*surface.VectorT.X + vertices[vertex[i]].Y*surface.VectorT.Y + vertices[vertex[i]].Z*surface.VectorT.Z + surface.DistT

				uv.u[1] = vertices[vertex[i+1]].X*surface.VectorS.X + vertices[vertex[i+1]].Y*surface.VectorS.Y + vertices[vertex[i+1]].Z*surface.VectorS.Z + surface.DistS
				uv.v[1] = vertices[vertex[i+1]].X*surface.VectorT.X + vertices[vertex[i+1]].Y*surface.VectorT.Y + vertices[vertex[i+1]].Z*surface.VectorT.Z + surface.DistT

				uv.u[0] = uv.u[0] / float32(mip.Width)
				uv.v[0] = 1.0 - uv.v[0]/float32(mip.Height)

				uv.u[1] = uv.u[1] / float32(mip.Width)
				uv.v[1] = 1.0 - uv.v[1]/float32(mip.Height)

				uv.u[2] = uv.u[2] / float32(mip.Width)
				uv.v[2] = 1.0 - uv.v[2]/float32(mip.Height)

				uv_list = append(uv_list, uv)

				for i := range uv.u {
					u_tmp := math.Abs(float64(uv.u[i]))
					v_tmp := math.Abs(float64(uv.v[i]))

					if u_tmp < u_min {
						u_min = u_tmp
					}
					if u_tmp > u_max {
						u_max = u_tmp
					}
					if v_tmp < v_min {
						v_min = v_tmp
					}
					if v_tmp > v_max {
						v_max = v_tmp
					}
				}

				// FIXME: calculate correct UV
				uv_light_list = append(uv_light_list, 0)
				uv_light_list = append(uv_light_list, 1)
				uv_light_list = append(uv_light_list, 0)
				uv_light_list = append(uv_light_list, 1)
				uv_light_list = append(uv_light_list, 1)
				uv_light_list = append(uv_light_list, 1)
			}

			if face.Lightmap != -1 {

				var light_w = int(((u_max-u_min) * float64(mip.Width))/16.0+2.0)
				var light_h = int(((v_max-v_min) * float64(mip.Height))/16.0+2.0)

		//		light_w = int((v_max-v_min) * float64(mip.Width)/16.0)
		//		light_h = int((u_max-u_min) * float64(mip.Height)/16.0)

				fmt.Printf("--> %d %d\n", mip.Width, mip.Height)
				fmt.Printf("->u: %f, %f\n", u_max-u_min, (u_max-u_min) * float64(mip.Width))
				fmt.Printf("->v: %f, %f\n", v_max-v_min, (v_max-v_min) * float64(mip.Height))


				fmt.Printf("-> %dx%d\n", light_w, light_h)

				myimage := image.NewGray(image.Rect(0, 0, light_w, light_h))

				for i := 0; i < light_h; i++ {
					for j := 0; j < light_w; j++ {
						myimage.Set(j, i, color.Gray{lightmaps[face.Lightmap+int32(i*light_w+j)]})
					}
				}

				lightmap_images = append(lightmap_images, myimage)

				myfile, _ := os.Create("test.png")

				png.Encode(myfile, myimage)

				if (light_w > 10 && light_h > 10) {
					return
				}
			}
		}
		return
		var out *os.File

		if model_idx == 0 {
			out, err = os.Create(os.Args[2])
		} else {
			out, err = os.Create(strings.TrimSuffix(os.Args[2], filepath.Ext(os.Args[2])) + "_" + strconv.Itoa(model_idx) + ".json")
		}
		if err != nil {
			log.Fatal(err)
		}
		defer out.Close()

		out.WriteString("{\n")
		out.WriteString("\t\"materials\": [")

		for _, material := range material_list {
			out.WriteString("\n\t{\n")
			out.WriteString("\t\t\"mapDiffuse\": \"textures/" + strings.TrimRight(string(miptex[material].Name[:]), "\x00") + ".jpg\",\n")
			out.WriteString("\t\t\"mapDiffuseWrap\": [\"repeat\", \"repeat\"]\n")
			// out.WriteString("\t\t\"mapLight\": \"lightmap.png\"\n")
			out.WriteString("\t},")
		}
		out.Seek(-1, io.SeekCurrent)

		out.WriteString(" ],\n")
		out.WriteString("\t\"vertices\": [")

		for _, vertice := range vertice_list {
			out.WriteString(strconv.FormatFloat(float64(vertices[vertice].X), 'f', -1, 32) + ",")
			out.WriteString(strconv.FormatFloat(float64(vertices[vertice].Y), 'f', -1, 32) + ",")
			out.WriteString(strconv.FormatFloat(float64(vertices[vertice].Z), 'f', -1, 32) + ",")
		}
		out.Seek(-1, io.SeekCurrent)

		out.WriteString("],\n")
		out.WriteString("\t\"uvs\": [[")

		for _, uv := range uv_list {
			out.WriteString(strconv.FormatFloat(float64(uv.u[0]), 'f', -1, 32) + ",")
			out.WriteString(strconv.FormatFloat(float64(uv.v[0]), 'f', -1, 32) + ",")
			out.WriteString(strconv.FormatFloat(float64(uv.u[1]), 'f', -1, 32) + ",")
			out.WriteString(strconv.FormatFloat(float64(uv.v[1]), 'f', -1, 32) + ",")
			out.WriteString(strconv.FormatFloat(float64(uv.u[2]), 'f', -1, 32) + ",")
			out.WriteString(strconv.FormatFloat(float64(uv.v[2]), 'f', -1, 32) + ",")
		}
		out.Seek(-1, io.SeekCurrent)

		out.WriteString("],[")

		// for _, uv := range uv_light_list {
		// 	out.WriteString(strconv.FormatFloat(float64(uv), 'f', -1, 32) + ",")
		// }
		// out.Seek(-1, io.SeekCurrent)

		out.WriteString("]],\n")
		out.WriteString("\t\"faces\": [")

		for i, face := range face_list {
			out.WriteString("10,")
			out.WriteString(strconv.Itoa(face.vert[0]) + ",")
			out.WriteString(strconv.Itoa(face.vert[1]) + ",")
			out.WriteString(strconv.Itoa(face.vert[2]) + ",")
			out.WriteString(strconv.Itoa(face.mat) + ",")
			out.WriteString(strconv.Itoa(3*i) + ",")
			out.WriteString(strconv.Itoa(3*i+1) + ",")
			out.WriteString(strconv.Itoa(3*i+2) + ",")
			// out.WriteString(strconv.Itoa(3*i) + ",")
			// out.WriteString(strconv.Itoa(3*i+1) + ",")
			// out.WriteString(strconv.Itoa(3*i+2) + ",")
		}
		out.Seek(-1, io.SeekCurrent)

		out.WriteString("]\n")
		out.WriteString("}\n")
	}
}
