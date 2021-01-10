from django import forms
from django.core.validators import MaxValueValidator
from django.forms import MultiWidget, NumberInput


class LatitudeWidget(forms.MultiWidget):
    def __init__(self, attrs=None):
        super().__init__([
            forms.NumberInput( attrs={'max':90, 'min':0, 'class': 'numinp'}),
            forms.NumberInput(attrs={'max':59, 'min':0, 'class': 'numinp'}),
            forms.NumberInput(attrs={'max':59, 'min':0, 'class': 'numinp'}),
            forms.Select(choices=(
            ('north', 'N'),
            ('south', 'S')), attrs={'class':'selectdrop'})])
        
    
    def decompress(self, value):
        if value:
            return value
        return [0, 0, 0, 'N']

class LatitudeField(forms.MultiValueField):
    widget = LatitudeWidget

    def __init__(self, *args, **kwargs):
        fields = (
            forms.IntegerField(label="lat0" ),
            forms.IntegerField(),
            forms.IntegerField(),
            forms.ChoiceField(choices=[('north', 'N'),('south', 'S')])
        )

        super().__init__(fields, *args, **kwargs)
    
    def compress(self, data_list):
        if data_list[3] == 'north':
            return data_list[0]+data_list[1]/60+data_list[2]/3600
        return -(data_list[0]+data_list[1]/60+data_list[2]/3600)

class LongitudeWidget(forms.MultiWidget):
    def __init__(self, attrs=None):
        super().__init__([
            forms.NumberInput(attrs={'max':180, 'min':0, 'class': 'numinp'}),
            forms.NumberInput(attrs={'max':59, 'min':0, 'class': 'numinp'}),
            forms.NumberInput(attrs={'max':59, 'min':0, 'class': 'numinp'}),
            forms.Select(choices=(
            ('est', 'E'),
            ('west', 'O')), attrs={'class':'selectdrop'})])
        
    
    def decompress(self, value):
        if value:
            return value
        return [0, 0, 0, 'E']

class LongitudeField(forms.MultiValueField):
    widget = LongitudeWidget

    def __init__(self, *args, **kwargs):
        fields = (
            forms.IntegerField(),
            forms.IntegerField(),
            forms.IntegerField(),
            forms.ChoiceField(choices=[('est', 'E'),('west', 'O')])
        )

        super().__init__(fields, *args, **kwargs)
    
    def compress(self, data_list):
        if data_list[3] == 'est':
            return data_list[0]+data_list[1]/60+data_list[2]/3600
        return -(data_list[0]+data_list[1]/60+data_list[2]/3600)

class parametersWidget(forms.MultiWidget):
    def __init__(self, attrs=None):
        super().__init__([
            forms.NumberInput(),
            forms.NumberInput()])
        
    
    def decompress(self, value):
        if value:
            return value
        return [0, 0]

class parametersField(forms.MultiValueField):
    widget = parametersWidget

    def __init__(self, *args, **kwargs):
        fields = ([
            forms.FloatField(),
            forms.FloatField()])
        

        super().__init__(fields, *args, **kwargs)
    
    def compress(self, data_list):
        if data_list:
            return data_list
        return [None, None]
class finalform(forms.Form):
    ellipsoid = forms.ChoiceField(choices=[('wgs', 'WGS84'),('grs', 'GRS80'), ('clarke', 'Clarke1880')], required=False, widget= forms.Select(attrs={'class':'selectdrop'}))
    grand =  forms.FloatField(required=False, min_value=0, initial=0, widget=forms.NumberInput(attrs={'class':'ccst'}))
    petit = forms.FloatField(required=False, min_value=0, initial=0, widget=forms.NumberInput(attrs={'class':'ccst'}))
    latitude = LatitudeField()
    longitude = LongitudeField()
    azimut = forms.FloatField(initial=0, widget=forms.NumberInput(attrs={'class':'ccst'}))
    distance_geodesique = forms.FloatField(min_value=0, initial=0, widget=forms.NumberInput(attrs={'class':'ccst'}))
    