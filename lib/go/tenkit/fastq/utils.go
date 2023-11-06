package fastq

import (
	"container/list"
)

type OrderedMapElement struct {
	value interface{}
	key   string
}

type OrderedMap struct {
	m map[string]*list.Element
	l list.List
}

func (m *OrderedMap) Init(size int) {
	m.l.Init()
	m.m = make(map[string]*list.Element, size)
}

func (self *OrderedMap) HasKey(key string) bool {
	_, ok := self.m[key]
	return ok
}

func (self *OrderedMap) Get(key string) (interface{}, bool) {
	if elem := self.m[key]; elem != nil {
		return elem.Value.(OrderedMapElement).value, true
	}
	return nil, false
}

func (self *OrderedMap) Pop() (string, interface{}, bool) {
	back := self.l.Back()
	if back == nil {
		return "", nil, false
	}

	elem := self.l.Remove(back).(OrderedMapElement)
	delete(self.m, elem.key)
	return elem.key, elem.value, true
}

func (self *OrderedMap) Set(key string, value interface{}) {
	if elem, ok := self.m[key]; ok {
		self.l.MoveToFront(elem)
	} else {
		elem := self.l.PushFront(
			OrderedMapElement{
				key:   key,
				value: value,
			})
		self.m[key] = elem
	}
}

func (self *OrderedMap) Length() int {
	return len(self.m)
}

func (self *OrderedMap) Values() []interface{} {
	values := make([]interface{}, 0, self.Length())
	for _, elem := range self.m {
		values = append(values, elem.Value.(OrderedMapElement).value)
	}
	return values
}
