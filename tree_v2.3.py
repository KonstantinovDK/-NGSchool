#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""
import sys
from ete3 import Tree, SeqMotifFace, TreeStyle, add_face_to_node, NodeStyle, TextFace, AttrFace, faces, NodeStyle 
from Bio import SeqIO
from Bio import Phylo
import os
import re


def atoi(text):
    return int(text) if text.isdigit() else text
def natural_keys(text):
    return [ atoi(c) for c in re.split('(\d+)', text) ]

class node_collor:
    #функция забора характеристик узла
    def __init__(self, name, color="",pfam=""):
        self.name   = name 
        self.color   = color
        self.pfam   = pfam
    #функция добавления color к вершине
    def new_color_(self, color):
        self.color=color
        #функция добавления color к вершине
    def new_pfam_(self, pfam):
        self.pfam=pfam
    #функция возвращающая описание узла
    def __str__(node):
        return (node.name+"\tcolor-"+node.color+"\tpfam-"+node.pfam+"\n")   
a=node_collor("node1")
print(a)
a.new_pfam_("thy")
print(a)
def layout(node):
    if node.is_leaf():
        N = AttrFace("name", fsize=14)
        faces.add_face_to_node(N, node, 1, position="branch-right")
#1__здесь можно уточнить клады
adres=os.getcwd()
all_clad=[]
file_Groups = open(adres+"/for_pic/Groups.txt", 'r')
for line in file_Groups:
    all_clad.append(line[:-1].split("\t"))
collor_list=["LightSteelBlue", "Moccasin","DarkSeaGreen","LemonChiffon" ,"#ffcccc" ,"Teal", "Aqua","Navy","Violet","Salmon","Pink"]
#2__здесь можно уточнить приорететы укоренения
ancestor_grop=[]
file_ = open(adres+"/for_pic/ancestor.txt", 'r')
for line in file_:
    ancestor_grop.append(line.replace("\n",""))

#2__здесь можно указать порядок цветов доменов
dic_domain_pic_pic=["Olive","brown","LightSteelBlue", "Moccasin","DarkSeaGreen","LemonChiffon" ,"#ffcccc" ,"Teal", "Aqua","orange", "green","pink","blue","yellow","red","Purple","Olive","brown","LightSteelBlue", "Moccasin","DarkSeaGreen","LemonChiffon" ,"#ffcccc" ,"Teal", "Aqua","orange", "green","pink","blue","yellow","red","Purple"]

#3_здесь можно указать какие узлы стоит удалить
dell_node=["PPE","OSI"]

#3_здесь можно указать какие узлы стоит сохранить
save_node=["AT","SL","OL","OS","ZM","PP"]


list_legend=["a1","a2","a3","a4","a5","a6","a7","a8","a9","a11","a12","a13","a14","a15","a16","a17","a18","a19","a20","a21","a22","a23",]
seq_seq="MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM"
files_all_in = os.listdir(adres+"/for_pic/1_tree_nwk/") 
for file in files_all_in:
    #print(file)
    if not(file.startswith(".")):
        #print(file)
        Phylo.convert(adres+"/for_pic/1_tree_nwk/"+file, 'newick', adres+"/for_pic/4_tree_nwk/"+file, 'newick')

def get_example_tree(File):
    adres=os.getcwd()
    file_out_supliment = open(adres+"/out_spliment/"+File, 'w')
    node_file = open(adres+"/node/"+File, 'w')
    # Create a random tree and add to each leaf a random set of motifs
    # from the original set
    #t = Tree("( (A, B, C, D, E, F, G), H, I);")
    #Считываем все домены
    domain_all_legend={}
    file_all_domen=os.listdir(adres+"/for_pic/1_tree_nwk/") 
    file_all_domen.remove(".DS_Store")
    file_all_domen.sort()
    i=0
    for file_domain in file_all_domen:
        file_open_domain = open(adres+"/for_pic/3_domain/"+file_domain, 'r')
        for line in file_open_domain:
            line_=line.split("\t")
            try:           
                if not (line_[2] in domain_all_legend):
                    domain_all_legend.setdefault(line_[2],dic_domain_pic_pic[i])
                    i+=1
                if i>len(dic_domain_pic_pic):
                    i=0
            except:
                a=0


    mem=""
    file_open = open(adres+"/for_pic/1_tree_nwk/"+File, 'r')
    for line in file_open:
        mem=mem+line
    tt = Tree(mem, format=0)   


    style = NodeStyle()   
    style["fgcolor"] = "#000000"
    style["size"] = 0
    style["vt_line_color"] = "#000000"
    style["hz_line_color"] = "#000000"
    style["vt_line_width"] = 4
    style["hz_line_width"] = 4
    style["vt_line_type"] = 8 # 0 solid, 1 dashed, 2 dotted
    style["hz_line_type"] = 8
    for node in tt.traverse("levelorder"):
        node.img_style = style
        if (len(node.name))>1:
            node_file.write(node.name+"\n")
        children1=node.children
        for element in children1:
            element.img_style = style
            
                
    for node in tt.traverse("preorder"):
        node.img_style = style
        children1=node.children
        for element in children1:
        	element.img_style = style
    node_file.close

    #вывести дерево с цветами
    #print (tt.get_ascii(attributes=["name", "color"], show_internal=False))
    #поиск предка
    ancestor1=""
    i=0
    for element in ancestor_grop:
        if i==0:
            for node in tt.traverse("postorder"):
                if i==0:
                    node_name=str(node.name)
                    if (node_name.startswith(element)) and  not((node_name.startswith("PPE"))):
                        ancestor1=str(node.name)
                        i=1
                        break
                else:
                    break
        else:
            break
    if not (ancestor1==""):
        tt.set_outgroup(ancestor1)
        #tt.render(adres+"/out/"+File[:-3]+"_2.png", tree_style=circular_style)
        print(str(ancestor1)+" - предок")
        file_out_supliment.write(str(ancestor1)+"\t"+" - предполагаемый корень"+"\n")
    else:
        print("Не нашел предка")



    file_out_supliment.write("\n\n\n Выявленные клады\n")
    #добавляем цвета к кладам
    for leaf in tt:
        i=0
        node_name=str(leaf.name)
        for clad in all_clad:
            collor=collor_list[i]
            i+=1
            for element in clad:
                if (node_name.startswith(element)):
                    leaf.add_features(color=collor)
                    #print(leaf)
    #print(tt)
    #забираем монофилитические цвета
    #print (tt.get_ascii(attributes=["name", "color"], show_internal=False))
    ii=-1
    for clad in all_clad:
        ii+=1
        collor=collor_list[ii]
        for monophyletic_tree in tt.get_monophyletic(values=[collor], target_attr="color"):
            i=[]
            name_node_mono_color=[]
            for leaf in monophyletic_tree:
                i.append(leaf)
                name_node_mono_color.append(leaf.name)
            if len(i)>1:
                n1 = tt.get_common_ancestor(i)
                nst1 = NodeStyle()
                nst1["bgcolor"] = collor
                nst1["fgcolor"] = "#000000"
                nst1["size"] = 0
                nst1["vt_line_color"] = "#000000"
                nst1["hz_line_color"] = "#000000"
                nst1["vt_line_width"] = 4
                nst1["hz_line_width"] = 4
                nst1["vt_line_type"] = 8 # 0 solid, 1 dashed, 2 dotted
                nst1["hz_line_type"] = 8
                n1.set_style(nst1)

                for element in name_node_mono_color:
                    file_out_supliment.write(str(element)+"\t"+" - "+collor+"\n")
                file_out_supliment.write("\n")

   
    file_out_supliment.write("\n\n\n Легенда доменного состава\n")
    #добавляем разметку по доменам
    dic_seq={}
    dic_domain={}
    dic_domain_pic={}
    i=0
    list_legend_domain3=[]
    for node in tt.traverse("postorder"):
        #длины белков
        fasta_sequences=SeqIO.parse(open(adres+"/for_pic/2_MSA/"+File), "fasta")
        for element in fasta_sequences:
            if str(element.id)==str(node.name):
                dic_seq.setdefault(str(node.name),str(element.seq))
        #доменный состав
        a=[]
        file_domain = open(adres+"/for_pic/3_domain/"+File, 'r')
        for line in file_domain:
            line_=line.split("\t")
            if line_[0]==str(node.name):
                if not (line_[2]  in list_legend_domain3):
                    list_legend_domain3.append(line_[2])



                if not (line_[2] in dic_domain_pic):
                    dic_domain_pic.setdefault(line_[2],dic_domain_pic_pic[i])
                    i+=1
                    #print(dic_domain_pic[line_[2]])
                    #print(i)

                    a1=[int(line_[3]),int(line_[4]), "()", None, 15, "black", domain_all_legend[line_[2]], "arial|9|black|"+line_[2]]
                    a.append(a1)
                    dic_domain.setdefault(str(node.name),a)
                    file_out_supliment.write(line_[2]+"\t"+domain_all_legend[line_[2]]+"\n")
                else:
                    a1=[int(line_[3]),int(line_[4]), "()", None, 15, "black", domain_all_legend[line_[2]], "arial|9|black|"+line_[2]]
                    a.append(a1)
                    dic_domain.setdefault(str(node.name),a)

    for element in dic_domain:
        #print(str(element)+" "+ str(dic_domain[element]))
        try:
            seqFace = SeqMotifFace(seq=dic_seq[element], motifs=dic_domain[element], seq_format="line")
            (tt & element).add_face(seqFace, 0, "aligned")
        except:
            seqFace = SeqMotifFace(seq=dic_seq[element],  seq_format="line", gapcolor="red")
            (tt & element).add_face(seqFace, 0, "aligned")
            print("except")

    #Рисуем легенду
    circular_style = TreeStyle()
    circular_style.show_leaf_name = False
    circular_style.show_branch_length = True
    circular_style.show_branch_support = True
    circular_style.scale = 75
    circular_style.tree_width = 50
    file_domain.close
    file_domain = open(adres+"/for_pic/3_domain/"+File, 'r')
    list_legend_domain={}
    list_legend_domain2=[]
    #считали список доменов
    i=0
    for line in file_domain:
        line_=line.split("\t")
        try:
            if  not(line_[2] in list_legend_domain2):
                #print(line_[2])
                list_legend_domain2.append(line_[2])
                list_legend_domain.setdefault("a"+str(i),line_[2])
                i+=1
        except:
            print("не понял что это за домен")
    i=0
    #считываем легенду доменов 
    file_domain_legend2={}
    file_domain_legend = open(adres+"/domain_legend.txt", 'r')
    for line in file_domain_legend:
        line_=line.split("\t")
        aaa=line_[1].replace(" ","_")
        aaa=aaa.replace("(","_")
        aaa=aaa.replace(")","_")
        aaa=aaa.replace(",","_")
        aaa=aaa.replace(":","_")
        aaa=aaa.replace(".","_")
        file_domain_legend2.setdefault(line_[0],aaa.replace("\n",""))
    #N = AttrFace("name", fsize=12)
    #faces.add_face_to_node(N, node, 1, position="branch-right")


    #рисуем домены
    ww=""
    for element in file_domain_legend2:
        ww=ww+","+file_domain_legend2[element]
    ww="("+ww[1:]+");"
    tree_domen_all=Tree(ww)
    for element in file_domain_legend2:
        try:
            element2=domain_all_legend[element]
            a1=[10,90, "()", None, 15, "black", domain_all_legend[element], "arial|9|black|"+element]
            i+=1
            a=[]
            a.append(a1)
            seqFace = SeqMotifFace(seq=seq_seq, motifs=a, seq_format="line")
            #node_node="a"+str(i)
            node_node=file_domain_legend2[element]
            try:
                (tree_domen_all & node_node).add_face(seqFace, 0, "aligned")
            except:
                q=1
                print("не нашел узел")
        except:
                q=1
    circular_style.layout_fn = layout
    tree_domen_all.render(adres+"/out_legend_all.png", tree_style=circular_style)


    file_domain_out = open(adres+"/123123123.txt", 'w')
    w=""
    for element in list_legend_domain3:
        w=w+","+file_domain_legend2[element]
    w="("+w[1:]+");"
    tree_domen=Tree(w)

    for element in list_legend_domain3:
        file_domain_out.write(element+"\n")
        a1=[10,90, "()", None, 15, "black", domain_all_legend[element], "arial|9|black|"+element]
        i+=1
        a=[]
        a.append(a1)
        try:
            seqFace = SeqMotifFace(seq=seq_seq, motifs=a, seq_format="line")
            #node_node="a"+str(i)
            node_node=file_domain_legend2[element]
            (tree_domen & node_node).add_face(seqFace, 0, "aligned")
        except:
            #print("Закончились узлы легенды")
            k=0
    circular_style.layout_fn = layout
    tree_domen.render(adres+"/out_legend/"+File[:-4]+".png", tree_style=circular_style)


    #удаленние части узлов 
    for node in tt.traverse("postorder"):
        try:                
            seqFace = SeqMotifFace(seq=dic_seq[str(node.name)], motifs=dic_domain[str(node.name)], seq_format="line")
            (tt & node.name).add_face(seqFace, 0, "aligned")
            a=0
            if len(node.name)<2:
                a=1
            for element_save in save_node:
                if (node.name).startswith(element_save):
                    a=1
            for element_dell in dell_node:
                if (node.name).startswith(element_dell):
                    a=0
            if a==0:
                node.delete()
        except:
            if len(node.name)>0:
                seqFace = SeqMotifFace(seq=dic_seq[str(node.name)],  seq_format="line", gapcolor="red")
                (tt & node.name).add_face(seqFace, 0, "aligned")
                node.delete()
                d0=0
    #удаленние части узлов ЗАВЕРШЕНО
    #особые точки
    node_color=[]
    file_node_color = open(adres+"/for_pic/4_color_node/out_list_gene2.txt", 'r')
    for line in file_node_color:
        node_color.append(line.replace("\n",""))
    for node in tt.traverse("postorder"):
        if node.name in node_color:
            style = NodeStyle()   
            style["fgcolor"] = "Red"
            style["size"] = 9
            style["vt_line_color"] = "#000000"
            style["hz_line_color"] = "#000000"
            style["vt_line_width"] = 4
            style["hz_line_width"] = 4
            style["vt_line_type"] = 8 # 0 solid, 1 dashed, 2 dotted
            style["hz_line_type"] = 8
            node.set_style(style)
    file_out_supliment.close
    #забираем монофилитические цвета
    #print (tt.get_ascii(attributes=["name", "color"], show_internal=False))
    ii=-1
    for clad in all_clad:
        ii+=1
        collor=collor_list[ii]
        for monophyletic_tree in tt.get_monophyletic(values=[collor], target_attr="color"):
            i=[]
            name_node_mono_color=[]
            for leaf in monophyletic_tree:
                i.append(leaf)
                name_node_mono_color.append(leaf.name)
            if len(i)>1:
                n1 = tt.get_common_ancestor(i)
                nst1 = NodeStyle()
                nst1["bgcolor"] = collor
                nst1["fgcolor"] = "#000000"
                nst1["size"] = 0
                nst1["vt_line_color"] = "#000000"
                nst1["hz_line_color"] = "#000000"
                nst1["vt_line_width"] = 4
                nst1["hz_line_width"] = 4
                nst1["vt_line_type"] = 8 # 0 solid, 1 dashed, 2 dotted
                nst1["hz_line_type"] = 8
                n1.set_style(nst1)

                for element in name_node_mono_color:
                    file_out_supliment.write(str(element)+"\t"+" - "+collor+"\n")
                file_out_supliment.write("\n")

    return tt




if __name__ == '__main__':
    adres=os.getcwd()
    file_all=os.listdir(adres+"/for_pic/1_tree_nwk/") 
    for File in file_all:
        if not(File.startswith(".")):
            print(File)
            t = get_example_tree(File)
            #выбор стиля 
            circular_style = TreeStyle()
            circular_style.show_leaf_name = False
            circular_style.show_branch_length = True
            circular_style.show_branch_support = True
            circular_style.scale = 75
            circular_style.tree_width = 50
            #circular_style.rotation = 90
            #circular_style.extra_branch_line_type=(0) 
            #circular_style.guiding_lines_type= (0)
            #circular_style.title.add_face(TextFace(File, fsize=25), column=0)
            circular_style.layout_fn = layout
            t.render(adres+"/out/"+File[:-4]+".png", tree_style=circular_style)
            circular_style.mode = "r" # draw tree in circular mode
            circular_style.extra_branch_line_type=(2) 
            circular_style.guiding_lines_type= (2)
            circular_style.guiding_lines_color =("red")
            for n in t.traverse():
			   nstyle = NodeStyle()
			   nstyle["fgcolor"] = "red"
			   nstyle["size"] = 15
			   n.set_style(nstyle) 
			   #N = AttrFace("name", fsize=30)
               #faces.add_face_to_node(N, node, 0, position="aligned")
            circular_style.legend.add_face =(TextFace("0.5 support"), column=1)
            circular_style.legend.add_face =(CircleFace(10, "red"), column=0)
            circular_style.layout_fn = layout
            circular_style.legend.add_face(TextFace("0.5 support"), column=1)
            circular_style.title.add_face(TextFace(File, fsize=40), column=0)
            TextFace_=TextFace(text, ftype='Verdana', fsize=10, fgcolor='black', penwidth=0, fstyle='normal', tight_text=False, bold=False)
            t.render(adres+"/out/"+File[:-4]+"_c.png", tree_style=circular_style)
