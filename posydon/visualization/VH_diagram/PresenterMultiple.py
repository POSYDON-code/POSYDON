"""Print multiple binaries side by side horizontally or 'hierarchically'"""

__authors__ = [
    "Wène Kouarfate <Wene.Kouarfate@etu.unige.ch>"
]

from posydon.visualization.VHdiagram import DisplayMode
from posydon.visualization.VH_diagram.Presenter import Presenter, PresenterMode, get_max_distance,CaseInfos,equal_with_epsilon,get_event_state_filename,file_exist,get_star_state_filename
from posydon.visualization.VH_diagram.GraphVisualizer import columnTYPE,StateInfos,ConnectedItem, GraphVisualizerCase

from PyQt5.QtWidgets import QApplication
from PyQt5.QtCore import QTimer
import matplotlib.pyplot as plt
from IPython.display import Image, display
from enum import Enum, auto
import os

class VHdiagramm_m:
    """Handle a multiple or hierarchical VH diagram plot."""

    def __init__(
        self,
        filename,
        path="./",
        index=[0],
        *,
        count_dict=None,
        hierarchy=False,
        presentMode=PresenterMode.DIAGRAM,
        displayMode=DisplayMode.INLINE_B,
        figsize=(10, 8)
    ):
        """Initialize a VHdiagram_m instance."""
        self._app = (
            QApplication.instance()
        )
        if not self._app:
            self._app = QApplication([])
    
        #PRESENTING HIERARCHICALLY
        if hierarchy :
            self._presenter = Presenter_h(filename=filename, path=path, count_dict=count_dict)
            self._presenter.present_h(index, presentMode)
        #PRESENTING MULTIPLES SIDE BY SIDE
        else:
            self._presenter = Presenter_m(filename=filename, path=path, count_dict=count_dict)
            self._presenter.present_m(index, presentMode)

        if displayMode == DisplayMode.INLINE_B:
            #void	singleShot(int msec, Functor functor)
            #This static function calls functor after a given time interval.
            QTimer.singleShot(0, lambda: self._display_inline_b())
        elif displayMode == DisplayMode.INLINE_S:
            QTimer.singleShot(0, lambda: self._display_inline_s(figsize))

        self._app.exec_()
    
    #Return indexes witch have been sorted hierarchycally for tree view
    def get_sorted_index(self):
        if type(self._presenter) == Presenter_h :
            return [dataf.index[0] for dataf in self._presenter._columns_data]
        else :#return empty if Presenter_h has no been called
            return []

    def _display_inline_b(self):
        filepath = self._presenter.screen()
        self._presenter.close()
        display(Image(filename=filepath))

    def _display_inline_s(self, figsize):
        filepath = self._presenter.screen()
        self._presenter.close()
        image = plt.imread(filepath)

        plt.figure(figsize=figsize)
        #Display data as an image, i.e., on a 2D regular raster.
        plt.imshow(image, aspect="auto")
        plt.axis("off")
    
class Presenter_m(Presenter):
    """Handle multiple Diagram view"""
    def __init__(self,filename, path="./", count_dict=None):
        super(Presenter_m, self).__init__(filename, path)
        self._count_dict = count_dict
        
    def present_m(self, indexes, mode=PresenterMode.DETAILED):
        """Preset the binary."""
        self._present_mode = mode
        self._visualizer.options().set_showed_mode(self._present_mode)
        self._visualizer().reset()
        
        if self._count_dict== None :
            self._count_dict = {i: 1/len(indexes) for i in indexes}
    
        for i in indexes :
            self._prepare_corresponding_data(i)
            self._update_visualisation_m()#redefined

        self._main_window.show()

    def _update_visualisation_m(self):
        """Update the display according to current_index and present_mode."""
        if self._current_index is None:
            return

        if self._present_mode == PresenterMode.DETAILED:
            self._prepare_basic_columns_m()
            self._detailed_presentation(self._current_data)
        elif self._present_mode == PresenterMode.REDUCED:
            self._prepare_basic_columns_m()
            self._reduced_presentation(self._current_data)
        elif self._present_mode == PresenterMode.SIMPLIFIED:
            self._prepare_basic_columns_m()
            self._simplified_presentation(self._current_data)
        elif self._present_mode == PresenterMode.DIAGRAM:
            self._prepare_diagram_columns_m()
            self._digram_presentation(self._current_data)
    
    def _prepare_basic_columns_m(self):#THE REST USES 3 COLUMNS COL_SPAN=1
        """Add columns for Detailled/Reduced/Simplified View."""
        self._time_id = self._visualizer().add_column(columnTYPE.TIMELINE)
        self._visualizer().get_column(self._time_id).set_title(f"{self._current_index} : TIME")

        self._S1_id = self._visualizer().add_column(columnTYPE.CONNECTED)
        self._visualizer().get_column(self._S1_id).set_title("S1")

        self._event_id = self._visualizer().add_column(columnTYPE.CONNECTED)
        self._visualizer().get_column(self._event_id).set_title("EVENT/STATE")

        self._S2_id = self._visualizer().add_column(columnTYPE.CONNECTED)
        self._visualizer().get_column(self._S2_id).set_title("S2")
    
    def _prepare_diagram_columns_m(self):
        """Add columns for multiple Diagram View."""
        self._time_id = self._visualizer().add_column(columnTYPE.TIMELINE)
        self._state_id = self._visualizer().add_column(columnTYPE.CONNECTED)
        self._visualizer().get_column(self._state_id).set_title(f"{self._current_index} : {self._count_dict[self._current_index] : .2f} %")

class Presenter_h(Presenter):
    """Handles drawing binaries as hierarchical trees"""
    def __init__(self,filename, path="./", count_dict=None):
        #INITIALIZING PARENT
        super(Presenter_h, self).__init__(filename, path)
        self._count_dict = count_dict
        self._column_width = 400 #fixed but could be dynamic
        
    def present_h(self, indexes, mode=None, count_dict=None):
        """Preset the binary."""
        self._present_mode = PresenterMode.DIAGRAM
        self._visualizer.options().set_showed_mode(self._present_mode) #keep
        self._visualizer().reset()
        
        if self._count_dict == None :
            self._count_dict = {i: 1/len(indexes) for i in indexes}
        
        self._columns_indexes = []
        self._columns_data = []
        self._infos = []
        self._edges = []
    
        for i in indexes :
            self._prepare_corresponding_data(i)
            self._columns_data.append(self._current_data.copy())
            self._update_visualisation_h()
            
        self._diagram_presentation_h()
        self._main_window.show()
        
    def _sort_h(self):
        """Hierarchycally sorting binaries using python native hierachical
        sorting of variable length tuples provided by _tuple_h"""
        self._infos.sort(key=lambda info_col : self._tuple_h(info_col))
        #Reajusting column_id witch may be different from when created in _prepare_diagram_line()
        new_col_data = []
        
        for i in range(len(self._infos)):
            new_col_data.append(
                self._columns_data[self._infos[i][0].column_id]
            ) 
            for info in self._infos[i]:                   
                info.column_id = i
        
        self._columns_data = new_col_data

    def _tuple_h(self, info_col):
        return [(info.S1_filename.split("\\")[-1] +
                info.S2_filename.split("\\")[-1] if (
                    info.event_filename == None
                )else info.event_filename.split("\\")[-1]) if (
            type(info)==StateInfos
        )else info.centered_text
                for info in info_col]
    
    def _compareInfos(self, info0, info1):
        """Conpare two Infos objects (ca be CaseInfos or StateInfos)"""
        if type(info0) == type(info1):
            if type(info0) == StateInfos :
                if info0.event_filename == None == info1.event_filename :
                    return (
                        info0.S1_filename.split("\\")[-1]==info1.S1_filename.split("\\")[-1]
                    ) and (
                        info0.S2_filename.split("\\")[-1]==info1.S2_filename.split("\\")[-1]
                    )
                elif info0.event_filename != info1.event_filename :
                    return False
                else : return info0.event_filename == info1.event_filename
            else :
                return info0.centered_text == info1.centered_text
        else : return False
    
    def _remove_duplicate(self,start_col,end_col,line):
        """Remove consecutive horizontally duplicated state to get the tree view"""
        if start_col > end_col :
            return
        if start_col == end_col :
            for info in self._infos[start_col][line:]:
                j = info.column_id
                rel = self._count_dict[self._columns_data[j].index[0]]
                if type(info)==StateInfos :
                    info.top_texts = [f"{rel: .4f} %"]
                else :
                    info.top_left_text = f"{rel: .4f} %"

        shortest_col = min(self._infos[start_col : end_col+1] ,
                       key=lambda info_col: len(info_col))
        if len(shortest_col) > line:
            info_col0 = self._infos[start_col]
            
            #Range of columns with same state at current line
            rg_column = [[start_col]]        
            for j in range(start_col+1, end_col+1):
                info_col1 = self._infos[j]
                if self._compareInfos(info_col0[line], info_col1[line]):
                    #replacing by empty state
                    info_col0[line] = CaseInfos(
                        info_col0[line].column_id,
                        centered_txt = " "
                    )
                    info_col0[line].connected = False
                else:
                    #recording a new range
                    rg_column[-1].append(j-1)
                    rg_column.append([j])
                info_col0 = info_col1
            rg_column[-1].append(end_col)
            self._edges.extend([[line-1, start_col + (end_col-start_col)//2,
                                 line, 
                                 c[0] + (c[1]-c[0])//2] for c in rg_column])
            
            self._center_column(rg_column, line)
            
            for start,stop in rg_column :
                self._remove_duplicate(start, stop, line+1)
        else :
            no_col = self._infos.index(shortest_col, start_col, end_col+1)
            self._remove_duplicate(start_col, no_col-1, line)
            self._remove_duplicate(no_col+1, end_col, line)
            for i in range(start_col,end_col+1):
                if i == no_col or len(self._infos[i])<=line:
                    continue
                else:
                    if type(self._infos[i][line]) == StateInfos:
                        self._edges.append([line-1, start_col + (end_col-start_col)//2,
                                             line, i])
            
        
    def _center_column(self, range_column, line):
        for k,j in range_column:
            if j - k >= 1 :
                info = self._infos[j][line]
                self._infos[j][line] = CaseInfos(
                        info.column_id,
                        centered_txt = " "
                    )
                self._infos[j][line].connected = False
                
                info.column_id = self._infos[k + (j-k)//2][line].column_id
                rel = sum(
                    [self._count_dict[self._columns_data[i].index[0]]
                    for i in range(k,j+1)])
                if type(info)==StateInfos :
                    info.top_texts = [f"{rel: .4f} %"]
                else :
                    info.top_left_text = f"{rel: .4f} %"
                self._infos[k + (j-k)//2][line] = info
            elif j==k :
                info = self._infos[j][line]
                rel = self._count_dict[self._columns_data[j].index[0]]
                if type(info)==StateInfos :
                    info.top_texts = [f"{rel: .4f} %"]
                else :
                    info.top_left_text = f"{rel: .4f} %"
            
        

    def _update_visualisation_h(self):
        """Update the display according to current_index and present_mode."""
        if self._current_index is None:
            return

        if self._present_mode == PresenterMode.DETAILED:
            self._prepare_basic_columns()
        elif self._present_mode == PresenterMode.REDUCED:
            self._prepare_basic_columns()
            self._reduced_presentation(self._current_data)
        elif self._present_mode == PresenterMode.SIMPLIFIED:
            self._prepare_basic_columns()
            self._simplified_presentation(self._current_data)
        elif self._present_mode == PresenterMode.DIAGRAM:
            self._prepare_diagram_columns_h()
            
    
    def _prepare_diagram_columns_h(self):
        """Add columns for Diagram View."""
        self._state_id = self._visualizer().add_column(columnTYPE.CONNECTED)
        self._visualizer()._layout.setColumnMinimumWidth(self._state_id, self._column_width)
        self._columns_indexes.append(self._state_id)
    
    def _diagram_presentation_h(self):
        
        simplified_datas = []
        for i in self._columns_indexes:
            simplified_datas.append(
                self._simplify_data(self._columns_data[i])
            )
        
        for i in self._columns_indexes:
            simplified_data = simplified_datas[i]
            self._state_id = i
            max_distance_data = max(simplified_data, key=get_max_distance)
            max_distance = get_max_distance(max_distance_data)

            self._visualizer()._columns[i]._row_index = 16
            
            self._infos.append([])
            for line_data in simplified_data:
                dline = self._prepare_diagram_line(line_data, max_distance)
                self._infos[-1].append(dline[0])
            
            if (
                simplified_data[-1]["state"].state_after == "disrupted"
                or simplified_data[-1]["state"].state_after == "merged"
            ):
                aditional_info = CaseInfos(self._state_id)
                aditional_info.centered_text = simplified_data[-1][
                    "state"].state_after
                self._infos[-1].append(aditional_info)
        
        self._sort_h()
        self._remove_duplicate(0,len(self._infos)-1,0)
        
        for i in self._columns_indexes:
            self._state_id = i
            self._visualizer()._columns[i]._row_index = 16
            binary = self._columns_data[i].index[0]
            self._visualizer().get_column(i).set_title(f"{binary} : {self._count_dict[binary]: .4f} %")
            for line_info in self._infos[i]:
                self._visualizer().add_line([line_info])
           
        self._paint_obl_edges()
    
    def _prepare_diagram_line(self, data, max_distance):
        
        state_info = StateInfos(self._state_id)
        if self._distance_representation and "separation" in data:
            state_info.distance = self._get_distance_representation(
                data["separation"], max_distance
            )
        else:
            state_info.distance = 1

        state_before_filename = os.path.join(
            self.PATH_TO_DRAWS,
            get_event_state_filename(
                data["S1_state"].state_before,
                data["state"].state_before,
                data["S2_state"].state_before,
                suffix=".png",
            ),
        )
        state_after_filename = os.path.join(
            self.PATH_TO_DRAWS,
            get_event_state_filename(
                data["S1_state"].state_before,
                data["state"].state_before,
                data["S2_state"].state_before,
                suffix=".png",
            ),
        )
        if (data["state"].state_after != "detached"
                and file_exist(state_after_filename)):
            state_info.event_filename = state_after_filename
        elif data["state"].state_before != "detached" and file_exist(
            state_before_filename
        ):
            state_info.event_filename = state_before_filename
        else:
            if data["S1_state"].state_after is not None:
                state_info.S1_filename = os.path.join(
                    self.PATH_TO_DRAWS,
                    get_star_state_filename(
                        data["S1_state"].state_after, suffix=".png"
                    ),
                )
            else:
                state_info.S1_filename = os.path.join(
                    self.PATH_TO_DRAWS,
                    get_star_state_filename(
                        data["S1_state"].state_before, suffix=".png"
                    ),
                )

            if data["S2_state"].state_after is not None:
                state_info.S2_filename = os.path.join(
                    self.PATH_TO_DRAWS,
                    get_star_state_filename(
                        data["S2_state"].state_after, suffix=".png"
                    ),
                )
            else:
                state_info.S2_filename = os.path.join(
                    self.PATH_TO_DRAWS,
                    get_star_state_filename(
                        data["S2_state"].state_before, suffix=".png"
                    ),
                )
        state_info.connected = False #Disabling default vertical edge
            
        return [state_info]
    
    def _paint_obl_edges(self):
        """Add ConnectedItem objects corresponding to edges"""
        for i1,j1,i2,j2 in self._edges[1:] :
            if type(self._visualizer()._columns[j1]._items[i1]) == GraphVisualizerCase : continue
            self._visualizer()._columns[j1]._connected_items.append(
                        ConnectedItem(
                            self._visualizer()._columns[j1]._items[i1],
                            self._visualizer()._columns[j2]._items[i2]))